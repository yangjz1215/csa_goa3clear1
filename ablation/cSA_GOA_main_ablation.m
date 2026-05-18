function [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive, best_scalar_solution, best_utility_solution, best_knee_solution] = cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, variant)
    if nargin < 13 || isempty(variant)
        variant = 'proposed';
    end

    params.variant = variant;
    params.RRH = RRH;
    params.RRH_type = RRH_type;
    params.UAV_type = UAV_type;
    if ~isfield(params, 'enable_migration_log')
        params.enable_migration_log = false;
    end
    if ~isfield(params, 'migration_stagnation_iters') || isempty(params.migration_stagnation_iters)
        params.migration_stagnation_iters = 20;
    end
    if ~isfield(params, 'enable_pv_interpolation')
        params.enable_pv_interpolation = true;
    end
    if ~isfield(params, 'pv_interpolation_interval') || isempty(params.pv_interpolation_interval)
        params.pv_interpolation_interval = 15;
    end
    if ~isfield(params, 'pv_interpolation_min_archive') || isempty(params.pv_interpolation_min_archive)
        params.pv_interpolation_min_archive = 10;
    end
    if ~isfield(params, 'mem_quota_m') || isempty(params.mem_quota_m)
        params.mem_quota_m = 2;
    end
    if ~isfield(params, 'pv_mix_logit_k') || isempty(params.pv_mix_logit_k)
        params.pv_mix_logit_k = -5;
    end
    if ~isfield(params, 'pv_mix_logit_c') || isempty(params.pv_mix_logit_c)
        params.pv_mix_logit_c = 0.38;
    end
    params = configureAblationVariant(params, variant);

    % [DEBUG] 运行时验证变体标志是否正确生效
    fprintf('[ABLATION_DEBUG] variant=%s | enable_phi_t=%d | enable_pv_interpolation=%d | enable_elite_migration=%d | enable_multi_subpop=%d\n', ...
        variant, params.enable_phi_t, params.enable_pv_interpolation, params.enable_elite_migration, params.enable_multi_subpop);

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
        % 公平性：无多子群时，单种群使用 3*K 候选解以保持总评估量一致
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
    stagnation_counter = zeros(n_subpops, 1);
    prev_fits = zeros(n_subpops, 1);

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

    cached_phi_t = 1.0;
    cached_bestUAV = zeros(0);

    for iter = 2:params.FES_max
        t = 1 - iter / params.FES_max;

        % 序参量缓存：bestUAV不变则复用；no_phi_t 变体跳过计算
        if params.enable_phi_t
            if size(cached_bestUAV, 1) ~= size(bestUAV, 1) || isempty(cached_bestUAV) || any(cached_bestUAV(:) ~= bestUAV(:))
                cached_phi_t = computePhasePhi(iter, params.FES_max, bestUAV, User, priorities, params, RRH);
                cached_bestUAV = bestUAV;
            end
            phi_t = cached_phi_t;
        else
            phi_t = 1.0;
        end
        if params.enable_phi_t
            pv_accept = 1 / (1 + exp(-(params.pv_mix_logit_k * (phi_t - params.pv_mix_logit_c))));
        else
            pv_accept = 1.0;  % 无 φ_t 门控时 PV 交换 100% 触发
        end

        % [DEBUG] 仅在迭代2打印一次，验证 phi_t/pv_accept/n_subpops 实际值
        if iter == 2
            fprintf('[ABLATION_DEBUG] iter=%d phi_t=%.4f pv_accept=%.4f n_subpops=%d K=%d\n', ...
                iter, phi_t, pv_accept, n_subpops, params.K);
        end

        for g = 1:n_subpops
            capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
        end

        for g = 1:n_subpops
            candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
            candidates = zeros(params.K, N_UAV, 2);

            random_leader = [];
            if params.enable_random_global_leader
                rl = sampleCandidates(subpops{g}, 1, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
                random_leader = squeeze(rl(1, :, :));
                if size(random_leader, 1) == 2 && size(random_leader, 2) == N_UAV
                    random_leader = random_leader';
                end
            end

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
                        if params.enable_phi_t
                            q_eff = max(0.05, min(0.95, params.subpop_params.q(g) * (1 - 0.50 * phi_t)));
                        else
                            q_eff = params.subpop_params.q(g);
                        end
                        if rand >= q_eff
                            pos = goaUShape(subpops{g}, mem_ref_pos, t, X_init, g);
                        else
                            pos = goaVShape(subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                        end
                    else
                        pos = X_init;
                    end

                    pos = max(Lb, min(Ub, pos));
                    candidates(i, uav_idx, :) = pos(:)';
                end
            end

            for i = 1:params.K
                cand_i = squeeze(candidates(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                    cand_i = cand_i';
                end

                if params.enable_pareto_leader && length(pareto_archive) >= 3
                    arch_utils = [pareto_archive.Utility];
                    arch_lats = [pareto_archive.Latency];
                    arch_nrgs = [pareto_archive.Energy];

                    [~, max_u_idx] = max(arch_utils);
                    leader_G1 = reshape(pareto_archive(max_u_idx).UAV_pos, N_UAV, 2);
                    [~, min_l_idx] = min(arch_lats);
                    leader_G2 = reshape(pareto_archive(min_l_idx).UAV_pos, N_UAV, 2);
                    [~, min_e_idx] = min(arch_nrgs);
                    leader_G3 = reshape(pareto_archive(min_e_idx).UAV_pos, N_UAV, 2);

                    if g == 1
                        global_leader = leader_G1;
                    elseif g == 2
                        global_leader = leader_G2;
                    else
                        global_leader = leader_G3;
                    end
                elseif params.enable_random_global_leader && ~isempty(random_leader)
                    global_leader = random_leader;
                else
                    global_leader = bestUAV;
                end

                if params.enable_goa_turn
                    if params.enable_phi_t
                        cap_eff = capturability_g(g) * (0.65 + 0.35 * (1 - phi_t));
                    else
                        cap_eff = capturability_g(g);
                    end
                    for uav_idx = 1:N_UAV
                        pos = goaTurn(cand_i(uav_idx, :), global_leader(uav_idx, :), cap_eff, t);
                        pos = max(Lb, min(Ub, pos));
                        cand_i(uav_idx, :) = pos(:)';
                    end
                end

                candidates(i, :, :) = reshape(cand_i, 1, N_UAV, 2);
            end

            if ~params.enable_multi_subpop
                g_eval = 1;
            else
                g_eval = g;
            end

            mem_matrix{g} = updateMemory(mem_matrix{g}, candidates, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g_eval, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
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
            N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

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

        if mod(iter, 3) == 0 || iter == params.FES_max
            do_pv_mix = params.enable_pv_interpolation && params.enable_multi_subpop && ...
                iter >= 20 && mod(iter, params.pv_interpolation_interval) == 0 && ...
                length(pareto_archive) >= params.pv_interpolation_min_archive && ...
                rand < pv_accept;
            if do_pv_mix
                mixed_candidates = pvInterpolationExchange(subpops, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU, User, priorities, params);
                for mc = 1:length(mixed_candidates)
                    cand_pos = mixed_candidates(mc).UAV_pos;
                    [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(cand_pos, User, priorities, params);
                    [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, cand_pos, cand_util, cand_lat, cand_nrg);
                end
            end

            for g = 1:n_subpops
                for i = 1:size(mem_matrix{g}, 1)
                    candidate = squeeze(mem_matrix{g}(i, :, :));
                    if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                        candidate = reshape(candidate, N_UAV, 2);
                    end
                    [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(candidate, User, priorities, params);
                    [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
                end
            end
        end

        if params.enable_multi_subpop && params.enable_elite_migration
            % 【速度优化】：不再重复评估所有解，只评估子种群中心点来判断是否停滞
            curr_subpop_fits = zeros(1, n_subpops);
            for g = 1:n_subpops
                center_pos = squeeze(local_mus(g, :, :));
                if size(center_pos, 1) == 1 && size(center_pos, 2) == N_UAV * 2
                    center_pos = reshape(center_pos, N_UAV, 2);
                end
                [curr_subpop_fits(g), ~, ~, ~] = calcFitness(center_pos, User, priorities, ...
                    E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                    N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
            end

            for g = 1:n_subpops
                if curr_subpop_fits(g) > prev_fits(g) + 1e-6
                    stagnation_counter(g) = 0;
                    prev_fits(g) = curr_subpop_fits(g);
                else
                    stagnation_counter(g) = stagnation_counter(g) + 1;
                end

                if iter > 20 && stagnation_counter(g) >= params.migration_stagnation_iters
                    if params.enable_migration_log
                        fprintf('[精英迁移] 迭代 %d: 子种群 G%d 连续 %d 代无改进，触发精英迁移\n', ...
                            iter, g, stagnation_counter(g));
                    end
                    mem_matrix{g} = migrateElite(mem_matrix, g, Ub, Lb, User, priorities, ...
                        E_remaining, params.E_max, params.k_move, params.subpop_params, ...
                        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
                    stagnation_counter(g) = 0;

                    subpop_fits_new = zeros(1, size(mem_matrix{g}, 1));
                    for i = 1:size(mem_matrix{g}, 1)
                        candidate = squeeze(mem_matrix{g}(i, :, :));
                        if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                            candidate = reshape(candidate, N_UAV, 2);
                        end
                        [subpop_fits_new(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                            E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                            N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
                    end
                    prev_fits(g) = max(subpop_fits_new);
                end
            end
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
    params.enable_pareto_leader = true;
    params.enable_multi_subpop = true;
    params.enable_goa_repulsion = true;
    params.enable_goa_turn = true;
    params.enable_elite_migration = true;
    params.enable_random_global_leader = false;
    params.enable_pv_interpolation = true;
    params.enable_phi_t = true;  % φ_t 统一调度（核心创新）

    switch variant
        case 'proposed'
            % 全部开启，无修改
        case 'no_phi_t'
            % 关闭 φ_t 统一调度：GOA 用固定 q(g)、goaTurn 用固定 capturability、PV 交换概率恒为 1
            params.enable_phi_t = false;
        case {'no_pv_interpolation', 'no_pv_exchange'}
            params.enable_pv_interpolation = false;
        case {'no_migration', 'no_elite_migration'}
            params.enable_elite_migration = false;
        case {'no_subpop', 'no_multi_subpop'}
            % 多子群关闭时，PV 交换和精英迁移自动失效
            params.enable_multi_subpop = false;
            params.enable_elite_migration = false;
            params.enable_pv_interpolation = false;
        % --- 向后兼容：保留旧变体名 ---
        case {'no_pareto', 'no_pareto_leader'}
            params.enable_pareto_leader = false;
        case 'no_goa'
            params.enable_goa_repulsion = false;
            params.enable_goa_turn = false;
        case 'no_goa_repulsion'
            params.enable_goa_repulsion = false;
        case 'no_goa_turn'
            params.enable_goa_turn = false;
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
    ideal = min(norm_objs, [], 1);
    [~, knee_idx] = min(sqrt(sum((norm_objs - ideal).^2, 2)));

    selected_uav = pareto_archive(knee_idx).UAV_pos;
    solution = buildSolutionStruct(selected_uav, fallback_scalar_fitness, User, priorities, params, 'archive_knee');
    solution.Utility = pareto_archive(knee_idx).Utility;
    solution.Latency = pareto_archive(knee_idx).Latency;
    solution.Energy = pareto_archive(knee_idx).Energy;
end

function new_pos = goaUShape(subpop, mem_ref_pos, t, X_init, g)
    A_g = (2 * rand - 1) * [0.6, 0.7, 0.5];
    new_pos = X_init(:)' + 3 * cos(2 * pi * rand) * t * mean(subpop.sigma(:)) + A_g(g) * (mem_ref_pos(:)' - X_init(:)');
end

function new_pos = goaVShape(subpop, mem_ref_pos, t, X_init, X_mean, g)
    x = 2 * pi * rand;
    V_x = (x < pi) * (-x / pi + 1) + (x >= pi) * (x / pi - 1);
    B_g = (2 * rand - 1) * [0.5, 0.6, 0.4];
    new_pos = X_init(:)' + 3 * V_x * t * mean(subpop.sigma(:)) + B_g(g) * (X_mean(:)' - X_init(:)');
end

function new_pos = goaTurn(pos, global_best_uav, cap, t)
    delta = cap * norm(pos - global_best_uav);
    direction = sign(global_best_uav - pos);
    theta = randn * 0.2;
    rot_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    direction = direction * rot_matrix;
    new_pos = pos(:)' + t * delta .* direction;
end
