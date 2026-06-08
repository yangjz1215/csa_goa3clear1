function [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive, best_scalar_solution, best_utility_solution, best_knee_solution] = cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, variant)
    if nargin < 12 || isempty(variant)
        variant = 'proposed';
    end

    params.variant = variant;
    params.RRH = RRH;
    params.RRH_type = RRH_type;
    params.UAV_type = UAV_type;
    if ~isfield(params, 'mem_quota_m') || isempty(params.mem_quota_m)
        params.mem_quota_m = 2;
    end
    params = configureAblationVariant(params, variant);

    fprintf('[ABLATION] variant=%s | multi_subpop=%d | goa_turn=%d | goa_repulsion=%d | pareto_leader=%d\n', ...
        variant, params.enable_multi_subpop, params.enable_goa_turn, params.enable_goa_repulsion, params.enable_pareto_leader);

    subpops = initSubpopulations(N_UAV, User, RRH, priorities, params.subpop_params, Ub, Lb, params.cover_radius, params.D_RU);

    if ~params.enable_multi_subpop
        subpops = {subpops{1}};
        params.G_weights = 1.0;
        params.subpop_params.sigma0 = mean(params.subpop_params.sigma0, 1);
        params.subpop_params.sigma_min = mean(params.subpop_params.sigma_min(:));
        if isfield(params.subpop_params, 'w_inertia')
            params.subpop_params.w_inertia = mean(params.subpop_params.w_inertia(:));
        end
        if isfield(params.subpop_params, 'c')
            params.subpop_params.c = mean(params.subpop_params.c(:));
        end
        params.subpop_params.q = mean(params.subpop_params.q(:));
        params.subpop_params.beta = mean(params.subpop_params.beta(:));
        params.K = params.K * 3;
        n_subpops = 1;
    else
        n_subpops = 3;
    end

    if isfield(params, 'K')
        params.K = round(params.K);
        if ~isfield(params, 'max_K_cap') || isempty(params.max_K_cap)
            params.max_K_cap = 200;
        end
        params.K = max(10, min(params.max_K_cap, params.K));
    end

    mem_matrix = cell(n_subpops, 1);
    for g = 1:n_subpops
        mem_matrix{g} = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
    end

    cg_curve = zeros(1, params.FES_max);
    energy_consumption = zeros(1, params.FES_max);
    E_remaining = params.E_max * ones(N_UAV, 1);
    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

    capturability_g = zeros(n_subpops, 1);
    for g = 1:n_subpops
        capturability_g(g) = calcCapturability(subpops{g}, 1, params.FES_max, g);
    end

    if ~isfield(params, 'RRH_radius')
        params.RRH_radius = 120;
    end
    if ~isfield(params, 'enable_bilevel')
        params.enable_bilevel = true;
    end
    if ~isfield(params, 'B_total')
        params.B_total = 10e6;
    end
    if ~isfield(params, 'F_total')
        params.F_total = 10e9;
    end
    if ~isfield(params, 'max_latency')
        params.max_latency = 1.0;
    end
    if ~isfield(params, 'kappa')
        params.kappa = 1e-27;
    end
    if ~isfield(params, 'P_tx')
        params.P_tx = 1;
    end
    if ~isfield(params, 'noise')
        params.noise = 1e-13;
    end

    [initial_fit, initial_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
        User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
    scalar_best_fit = initial_fit;
    bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
        E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

    [real_u, real_l, real_e] = calcMEC_Objectives(bestUAV, User, priorities, params);
    cg_curve(1) = real_u;
    energy_consumption(1) = initial_energy;

    fprintf('初始化完成：综合适应度=%.4f | 真实效用(优先级和)=%.1f | 时延=%.2fs | 能耗=%.1f J\n', ...
        scalar_best_fit, real_u, real_l, real_e);

    % 适应度缓存：避免 updateMemory 和 calcGlobalFitness 重复评估
    mem_fits_cache = cell(n_subpops, 1);
    mem_utils_cache = cell(n_subpops, 1);
    mem_lats_cache = cell(n_subpops, 1);
    mem_nrgs_cache = cell(n_subpops, 1);

    % 自适应权重初始化
    adaptive_weights = [0.70, 0.15, 0.15;
                        0.30, 0.50, 0.20;
                        0.20, 0.15, 0.65];
    if ~isfield(params, 'enable_adaptive_weight')
        params.enable_adaptive_weight = true;
    end
    params.test_weights = adaptive_weights;

    for iter = 2:params.FES_max
        t = 1 - iter / params.FES_max;

        for g = 1:n_subpops
            capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
        end

        for g = 1:n_subpops
            candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
            candidates = zeros(params.K, N_UAV, 2);

            current_mem_size = size(mem_matrix{g}, 1);
            if current_mem_size < params.K
                additional_candidates = sampleCandidates(subpops{g}, params.K - current_mem_size, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
                mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
            elseif current_mem_size > params.K
                mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
            end

            for i = 1:params.K
                cand_i = squeeze(candidates_init(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                    cand_i = cand_i';
                end
                X_mean_g = mean(cand_i, 1);

                for uav_idx = 1:N_UAV
                    X_init = cand_i(uav_idx, :);
                    mem_idx = min(i, size(mem_matrix{g}, 1));
                    mem_candidate = squeeze(mem_matrix{g}(mem_idx, :, :));
                    if size(mem_candidate, 1) == 2 && size(mem_candidate, 2) == N_UAV
                        mem_candidate = mem_candidate';
                    end
                    mem_ref_pos = mem_candidate(uav_idx, :);

                    if params.enable_goa_repulsion
                        q_eff = params.subpop_params.q(g);
                        if rand >= q_eff
                            pos = goaUShape(subpops{g}, mem_ref_pos, t, X_init, g);
                        else
                            pos = goaVShape(subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                        end
                    else
                        pos = X_init;
                    end

                    pos = projectToFeasiblePosition(pos, X_init, cand_i, uav_idx, RRH, N_RRH, N_UAV, params, Ub, Lb);
                    candidates(i, uav_idx, :) = pos(:)';
                end
            end

            if params.enable_pareto_leader && length(pareto_archive) >= 3
                arch_utils = [pareto_archive.Utility];
                arch_lats = [pareto_archive.Latency];
                arch_nrgs = [pareto_archive.Energy];

                top_pct = max(0.01, 0.20 - 0.19 * (iter / params.FES_max));
                n_archive = length(pareto_archive);
                top_n = max(3, round(n_archive * top_pct));

                % 拥挤度感知的Leader选择：从top-N中选择目标空间中最稀疏的解
                [~, sort_u_idx] = sort(arch_utils, 'descend');
                top_u_idx = sort_u_idx(1:top_n);
                top_u_objs = zeros(top_n, 3);
                for j = 1:top_n
                    top_u_objs(j, :) = [-arch_utils(top_u_idx(j)), arch_lats(top_u_idx(j)), arch_nrgs(top_u_idx(j))];
                end
                nn_dist_u = inf(top_n, 1);
                for j = 1:top_n
                    d = sqrt(sum((top_u_objs - repmat(top_u_objs(j,:), top_n, 1)).^2, 2));
                    d(j) = inf;
                    nn_dist_u(j) = min(d);
                end
                [~, sp_u] = max(nn_dist_u);
                leader_G1 = reshape(pareto_archive(top_u_idx(sp_u)).UAV_pos, N_UAV, 2);

                [~, sort_l_idx] = sort(arch_lats, 'ascend');
                top_l_idx = sort_l_idx(1:top_n);
                top_l_objs = zeros(top_n, 3);
                for j = 1:top_n
                    top_l_objs(j, :) = [-arch_utils(top_l_idx(j)), arch_lats(top_l_idx(j)), arch_nrgs(top_l_idx(j))];
                end
                nn_dist_l = inf(top_n, 1);
                for j = 1:top_n
                    d = sqrt(sum((top_l_objs - repmat(top_l_objs(j,:), top_n, 1)).^2, 2));
                    d(j) = inf;
                    nn_dist_l(j) = min(d);
                end
                [~, sp_l] = max(nn_dist_l);
                leader_G2 = reshape(pareto_archive(top_l_idx(sp_l)).UAV_pos, N_UAV, 2);

                [~, sort_e_idx] = sort(arch_nrgs, 'ascend');
                top_e_idx = sort_e_idx(1:top_n);
                top_e_objs = zeros(top_n, 3);
                for j = 1:top_n
                    top_e_objs(j, :) = [-arch_utils(top_e_idx(j)), arch_lats(top_e_idx(j)), arch_nrgs(top_e_idx(j))];
                end
                nn_dist_e = inf(top_n, 1);
                for j = 1:top_n
                    d = sqrt(sum((top_e_objs - repmat(top_e_objs(j,:), top_n, 1)).^2, 2));
                    d(j) = inf;
                    nn_dist_e(j) = min(d);
                end
                [~, sp_e] = max(nn_dist_e);
                leader_G3 = reshape(pareto_archive(top_e_idx(sp_e)).UAV_pos, N_UAV, 2);

                if g == 1
                    global_leader = leader_G1;
                elseif g == 2
                    global_leader = leader_G2;
                else
                    global_leader = leader_G3;
                end
            else
                global_leader = bestUAV;
            end

            for i = 1:params.K
                cand_i = squeeze(candidates(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                    cand_i = cand_i';
                end

                if params.enable_goa_turn
                    cap_eff = capturability_g(g);
                    for uav_idx = 1:N_UAV
                        pos = goaTurn(cand_i(uav_idx, :), global_leader(uav_idx, :), cap_eff, t);
                        pos = projectToFeasiblePosition(pos, cand_i(uav_idx, :), cand_i, uav_idx, RRH, N_RRH, N_UAV, params, Ub, Lb);
                        cand_i(uav_idx, :) = pos(:)';
                    end
                end

                candidate_pos = cand_i;
                if ~isFeasibleCandidate(candidate_pos, RRH, N_RRH, N_UAV, params)
                    cand_init_i = squeeze(candidates_init(i, :, :));
                    if size(cand_init_i, 1) == 2 && size(cand_init_i, 2) == N_UAV
                        cand_init_i = cand_init_i';
                    end
                    candidate_pos = cand_init_i;
                end
                g_eval = g;
            end

            [mem_matrix{g}, mem_fits_cache{g}, mem_utils_cache{g}, mem_lats_cache{g}, mem_nrgs_cache{g}] = ...
                updateMemory(mem_matrix{g}, candidates, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g_eval, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params, ...
                mem_fits_cache{g}, mem_utils_cache{g}, mem_lats_cache{g}, mem_nrgs_cache{g});
        end

        local_mus = zeros(n_subpops, N_UAV, 2);
        for g = 1:n_subpops
            for uav_idx = 1:N_UAV
                local_mus(g, uav_idx, :) = mean(squeeze(mem_matrix{g}(:, uav_idx, :)), 1);
            end
            subpops{g} = updateSubpopPV(subpops{g}, mem_matrix{g}, squeeze(local_mus(g, :, :)), params.subpop_params, g, iter, params.FES_max, N_UAV);
        end

        [curr_fit_best, curr_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
            User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params, mem_fits_cache);

        if curr_fit_best > scalar_best_fit
            scalar_best_fit = curr_fit_best;
            bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
                E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
        end

        [scalar_util, scalar_lat, scalar_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
        cg_curve(iter) = scalar_util;
        energy_consumption(iter) = curr_energy;

        if mod(iter, 50) == 0 || iter == params.FES_max
            fprintf('迭代 %d/%d, 综合得分: %.4f | 当前标量最优解真实效用: %.1f | 时延: %.2f s | 能耗: %.2f J | 档案解数量: %d\n', ...
                iter, params.FES_max, scalar_best_fit, scalar_util, scalar_lat, scalar_nrg, length(pareto_archive));
        end

        % Pareto存档更新：使用updateMemory缓存的目标值，避免重复调用calcMEC_Objectives
        % 性能优化：每5代更新一次Pareto存档（原来每代更新导致运行时间暴增4-5倍）
        if mod(iter, 5) == 0 || iter == params.FES_max
            for g = 1:n_subpops
                for i = 1:size(mem_matrix{g}, 1)
                    candidate = squeeze(mem_matrix{g}(i, :, :));
                    if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                        candidate = reshape(candidate, N_UAV, 2);
                    end
                    cand_util = mem_utils_cache{g}(i);
                    cand_lat = mem_lats_cache{g}(i);
                    cand_nrg = mem_nrgs_cache{g}(i);
                    [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
                end
            end
        end

        % 自适应权重旋转：每20代根据Pareto存档分布微调子种群权重方向
        if params.enable_adaptive_weight && mod(iter, 20) == 0 && iter > 20 && iter < 0.8 * params.FES_max && length(pareto_archive) >= 10
            arch_utils_aw = [pareto_archive.Utility];
            arch_lats_aw = [pareto_archive.Latency];
            arch_nrgs_aw = [pareto_archive.Energy];
            n_arch_aw = length(pareto_archive);

            range_u_aw = max(arch_utils_aw) - min(arch_utils_aw); if range_u_aw < 1e-9, range_u_aw = 1; end
            range_l_aw = max(arch_lats_aw) - min(arch_lats_aw); if range_l_aw < 1e-9, range_l_aw = 1; end
            range_e_aw = max(arch_nrgs_aw) - min(arch_nrgs_aw); if range_e_aw < 1e-9, range_e_aw = 1; end

            nu_aw = (arch_utils_aw - min(arch_utils_aw)) / range_u_aw;
            nl_aw = (arch_lats_aw - min(arch_lats_aw)) / range_l_aw;
            ne_aw = (arch_nrgs_aw - min(arch_nrgs_aw)) / range_e_aw;

            objs_aw = [nu_aw(:), nl_aw(:), ne_aw(:)];
            crowd_aw = zeros(n_arch_aw, 1);
            for ci = 1:n_arch_aw
                d_aw = sqrt(sum((objs_aw - repmat(objs_aw(ci,:), n_arch_aw, 1)).^2, 2));
                d_aw(ci) = inf;
                crowd_aw(ci) = min(d_aw);
            end

            [~, sp_aw] = max(crowd_aw);
            target_dir = [nu_aw(sp_aw), 1 - nl_aw(sp_aw), 1 - ne_aw(sp_aw)];
            target_dir = target_dir / sum(target_dir);

            min_d_aw = inf; closest_g_aw = 1;
            for gi = 1:3
                d_aw = norm(adaptive_weights(gi,:) - target_dir);
                if d_aw < min_d_aw
                    min_d_aw = d_aw;
                    closest_g_aw = gi;
                end
            end

            progress_aw = iter / params.FES_max;
            alpha_aw = 0.1 * (1 - progress_aw)^1.5;
            adaptive_weights(closest_g_aw, :) = (1 - alpha_aw) * adaptive_weights(closest_g_aw, :) + alpha_aw * target_dir;

            w_aw = adaptive_weights(closest_g_aw, :);
            w_aw = max(0.05, min(0.90, w_aw));
            adaptive_weights(closest_g_aw, :) = w_aw / sum(w_aw);

            params.test_weights = adaptive_weights;
        end
    end

    [final_util, final_lat, final_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
    [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, bestUAV, final_util, final_lat, final_nrg);

    best_scalar_solution = buildSolutionStruct(bestUAV, scalar_best_fit, User, priorities, params, 'scalar_best');
    best_utility_solution = selectMaxUtilitySolution(pareto_archive, bestUAV, scalar_best_fit, User, priorities, params);
    best_knee_solution = selectKneeSolution(pareto_archive, bestUAV, scalar_best_fit, User, priorities, params);
    best_fit = best_utility_solution.Utility;
end

function params = configureAblationVariant(params, variant)
    params.enable_multi_subpop = true;
    params.enable_goa_turn = true;
    params.enable_goa_repulsion = true;
    params.enable_pareto_leader = true;

    switch variant
        case 'proposed'
            % 全部开启
        case {'no_subpop', 'no_multi_subpop'}
            params.enable_multi_subpop = false;
            params.enable_goa_turn = false;
        case 'no_goa_turn'
            params.enable_goa_turn = false;
        case 'no_goa_repulsion'
            params.enable_goa_repulsion = false;
        case {'no_pareto_leader', 'no_pareto'}
            params.enable_pareto_leader = false;
        case 'no_adaptive_weight'
            params.enable_adaptive_weight = false;
        otherwise
            error('Unknown ablation variant: %s', variant);
    end
end

function solution = buildSolutionStruct(uav_pos, scalar_fitness, User, priorities, params, label)
    [util, lat, nrg, success_rate] = calcMEC_Objectives(uav_pos, User, priorities, params);
    solution = struct( ...
        'label', label, ...
        'UAV_pos', uav_pos, ...
        'ScalarFitness', scalar_fitness, ...
        'Utility', util, ...
        'Latency', lat, ...
        'Energy', nrg, ...
        'SuccessRate', success_rate);
end

function solution = selectMaxUtilitySolution(pareto_archive, fallback_uav, fallback_scalar_fitness, User, priorities, params)
    if ~isempty(pareto_archive)
        arch_util = [pareto_archive.Utility];
        [~, max_idx] = max(arch_util);
        selected_uav = pareto_archive(max_idx).UAV_pos;
        solution = buildSolutionStruct(selected_uav, fallback_scalar_fitness, User, priorities, params, 'archive_max_utility');
        solution.Utility = pareto_archive(max_idx).Utility;
        solution.Latency = pareto_archive(max_idx).Latency;
        solution.Energy = pareto_archive(max_idx).Energy;
    else
        solution = buildSolutionStruct(fallback_uav, fallback_scalar_fitness, User, priorities, params, 'archive_max_utility_fallback');
    end
end

function solution = selectKneeSolution(pareto_archive, fallback_uav, fallback_scalar_fitness, User, priorities, params)
    if isempty(pareto_archive)
        solution = buildSolutionStruct(fallback_uav, fallback_scalar_fitness, User, priorities, params, 'archive_knee_fallback');
        return;
    end

    U = [pareto_archive.Utility]';
    L = [pareto_archive.Latency]';
    E = [pareto_archive.Energy]';
    objs = [-U, L, E];
    mins = min(objs, [], 1);
    maxs = max(objs, [], 1);
    norm_objs = (objs - mins) ./ (maxs - mins + 1e-9);
    % 理想点：使用实际最小值而非理论最大值，避免Energy归一化被压缩
    ideal = [min(-U), min(L), min(E)];
    ideal_norm = (ideal - mins) ./ (maxs - mins + 1e-9);
    % 加权切比雪夫距离：max(w_j * |f_j - z_ideal|)，选择各目标均无短板的折中解
    w = [1/3, 1/3, 1/3];
    rho = 0.001;
    cheb_dist = zeros(size(norm_objs, 1), 1);
    for i = 1:size(norm_objs, 1)
        gap = w .* abs(norm_objs(i, :) - ideal_norm);
        cheb_dist(i) = max(gap) + rho * sum(gap);
    end
    [~, knee_idx] = min(cheb_dist);

    selected_uav = pareto_archive(knee_idx).UAV_pos;
    solution = buildSolutionStruct(selected_uav, fallback_scalar_fitness, User, priorities, params, 'archive_knee');
    solution.Utility = pareto_archive(knee_idx).Utility;
    solution.Latency = pareto_archive(knee_idx).Latency;
    solution.Energy = pareto_archive(knee_idx).Energy;
end

function new_pos = goaUShape(subpop, mem_ref_pos, t, X_init, g)
    r2 = 2 * pi * rand;
    a_coeffs = [0.6, 0.7, 0.5];
    A_g = (2 * rand - 1) * a_coeffs(g);
    sigma_mean = mean(subpop.sigma(:));
    new_pos = X_init(:)' + 3 * cos(r2) * t * sigma_mean + A_g * (mem_ref_pos(:)' - X_init(:)');
end

function new_pos = goaVShape(subpop, mem_ref_pos, t, X_init, X_mean, g)
    x = 2 * pi * rand;
    if x < pi
        V_x = -x / pi + 1;
    else
        V_x = x / pi - 1;
    end
    b_coeffs = [0.5, 0.6, 0.4];
    B_g = (2 * rand - 1) * b_coeffs(g);
    sigma_mean = mean(subpop.sigma(:));
    new_pos = X_init(:)' + 3 * V_x * t * sigma_mean + B_g * (X_mean(:)' - X_init(:)');
end

function new_pos = goaTurn(pos, global_best_uav, cap, t)
    delta = cap * norm(pos - global_best_uav);
    delta_vec = global_best_uav(:)' - pos(:)';
    dist = norm(delta_vec);
    if dist > 1e-10
        direction = delta_vec / dist;
    else
        direction = zeros(1, 2);
    end
    theta = randn * 0.2;
    rot_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    direction = direction * rot_matrix;
    new_pos = pos(:)' + t * delta .* direction;
end

function pos = projectToFeasiblePosition(pos, fallback_pos, candidate_positions, uav_idx, RRH, N_RRH, N_UAV, params, Ub, Lb)
    pos = max(Lb, min(Ub, pos));
    valid_pos = true;

    for rrh_idx = 1:N_RRH
        if norm(pos - RRH(rrh_idx, :)) < params.D_RU
            valid_pos = false;
            break;
        end
    end

    if valid_pos
        for other_uav = 1:N_UAV
            if other_uav ~= uav_idx
                other_pos = candidate_positions(other_uav, :);
                if norm(pos - other_pos) < params.D_UU
                    valid_pos = false;
                    break;
                end
            end
        end
    end

    if ~valid_pos
        pos = fallback_pos;
    end
end

function tf = isFeasibleCandidate(candidate_pos, RRH, N_RRH, N_UAV, params)
    tf = true;

    for rrh_idx = 1:N_RRH
        dists_rrh = sqrt(sum((candidate_pos - RRH(rrh_idx, :)).^2, 2));
        if any(dists_rrh < params.D_RU)
            tf = false;
            return;
        end
    end

    for uav_a = 1:N_UAV
        for uav_b = uav_a + 1:N_UAV
            if norm(candidate_pos(uav_a, :) - candidate_pos(uav_b, :)) < params.D_UU
                tf = false;
                return;
            end
        end
    end
end
