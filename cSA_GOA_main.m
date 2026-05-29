function [best_fit, bestUAV, cg_curve, energy_consumption, E_remaining_history, final_E_remaining, curr_curve, actual_iter, weighted_best_curve, pareto_archive, best_scalar_solution, best_utility_solution] = ...
    cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)

if ~isfield(params, 'enable_early_stop')
    params.enable_early_stop = true;
end
if ~isfield(params, 'enable_smart_stop')
    params.enable_smart_stop = true;
end
if ~isfield(params, 'enable_migration_log')
    params.enable_migration_log = false;
end
if ~isfield(params, 'enable_bilevel')
    params.enable_bilevel = true;
end
if ~isfield(params, 'enable_multi_subpop')
    params.enable_multi_subpop = true;
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
    params.pv_mix_logit_k = -5;  % 缓和logit门控曲线，保持中期PV信息交换频率
end
if ~isfield(params, 'pv_mix_logit_c') || isempty(params.pv_mix_logit_c)
    params.pv_mix_logit_c = 0.38;
end

params.RRH = RRH;
params.RRH_type = RRH_type;
params.UAV_type = UAV_type;

subpops = initSubpopulations(N_UAV, User, RRH, priorities, params.subpop_params, Ub, Lb, params.cover_radius, params.D_RU);

if isfield(params, 'K')
    params.K = max(10, round(params.K));
end

mem_matrix = cell(3, 1);
for g = 1:3
    mem_matrix{g} = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
end

cg_curve = zeros(1, params.FES_max);
curr_curve = zeros(1, params.FES_max);
energy_consumption = zeros(1, params.FES_max);
E_remaining_history = zeros(params.FES_max, N_UAV);
E_remaining = params.E_max * ones(N_UAV, 1);
weighted_best_curve = zeros(1, params.FES_max);
pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

stagnation_counter = zeros(3, 1);
prev_fits = zeros(3, 1);
capturability_g = zeros(3, 1);
for g = 1:3
    capturability_g(g) = calcCapturability(subpops{g}, 1, params.FES_max, g);
end

if ~isfield(params, 'RRH_radius')
    params.RRH_radius = 120;
end

[initial_fit, initial_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
    User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
    N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

scalar_best_fit = initial_fit;
bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
    E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

[~, weighted_initial, ~, ~] = calcFitness(bestUAV, User, priorities, ...
    E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
    N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
[init_util, init_lat, init_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);

cg_curve(1) = init_util;
curr_curve(1) = initial_fit;
weighted_best_curve(1) = weighted_initial;
energy_consumption(1) = initial_energy;
E_remaining_history(1, :) = E_remaining';

for g = 1:3
    subpop_fits = zeros(1, size(mem_matrix{g}, 1));
    for i = 1:size(mem_matrix{g}, 1)
        candidate = squeeze(mem_matrix{g}(i, :, :));
        if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
            candidate = reshape(candidate, N_UAV, 2);
        end
        [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
            E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
    end
    prev_fits(g) = max(subpop_fits);
end

fprintf('初始化完成：综合适应度=%.4f | 真实效用(优先级和)=%.1f | 时延=%.2fs | 能耗=%.1f J\n', ...
    scalar_best_fit, init_util, init_lat, init_nrg);

mo_stagnation_counter = 0;
actual_iter = params.FES_max;
cached_phi_t = 1.0;   % 序参量缓存：仅在bestUAV变化时重算
cached_bestUAV = zeros(0);

for iter = 2:params.FES_max
    t = 1 - iter / params.FES_max;

    % 序参量缓存：bestUAV不变则复用，避免每代重复调用calcMEC_Objectives
    if size(cached_bestUAV, 1) ~= size(bestUAV, 1) || isempty(cached_bestUAV) || any(cached_bestUAV(:) ~= bestUAV(:))
        cached_phi_t = computePhasePhi(iter, params.FES_max, bestUAV, User, priorities, params, RRH);
        cached_bestUAV = bestUAV;
    end
    phi_t = cached_phi_t;
    pv_accept = 1 / (1 + exp(-(params.pv_mix_logit_k * (phi_t - params.pv_mix_logit_c))));

    for g = 1:3
        capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
    end

    for g = 1:3
        candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
        candidates = zeros(params.K, N_UAV, 2);

        current_mem_size = size(mem_matrix{g}, 1);
        if current_mem_size < params.K
            additional_needed = params.K - current_mem_size;
            additional_candidates = sampleCandidates(subpops{g}, additional_needed, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
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

                q_eff = max(0.05, min(0.95, params.subpop_params.q(g) * (1 - 0.50 * phi_t)));  % 保留子群差异化分工
                if rand >= q_eff
                    pos = goaUShape(subpops{g}, mem_ref_pos, t, X_init, g);
                else
                    pos = goaVShape(subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                end

                pos = projectToFeasiblePosition(pos, X_init, cand_i, uav_idx, RRH, N_RRH, N_UAV, params, Ub, Lb);
                candidates(i, uav_idx, :) = pos(:)';
            end
        end

        if length(pareto_archive) >= 3
            arch_utils = [pareto_archive.Utility];
            arch_lats = [pareto_archive.Latency];
            arch_nrgs = [pareto_archive.Energy];

            % 动态 Top-% Leader 选择：前期从 Top 20% 随机选，后期收缩到 Top 1%
            top_pct = max(0.01, 0.20 - 0.19 * (iter / params.FES_max));
            n_archive = length(pareto_archive);
            top_n = max(3, round(n_archive * top_pct));

            % G1: Utility (越大越好)
            [~, sort_u_idx] = sort(arch_utils, 'descend');
            leader_G1 = reshape(pareto_archive(sort_u_idx(randi(top_n))).UAV_pos, N_UAV, 2);

            % G2: Latency (越小越好)
            [~, sort_l_idx] = sort(arch_lats, 'ascend');
            leader_G2 = reshape(pareto_archive(sort_l_idx(randi(top_n))).UAV_pos, N_UAV, 2);

            % G3: Energy (越小越好)
            [~, sort_e_idx] = sort(arch_nrgs, 'ascend');
            leader_G3 = reshape(pareto_archive(sort_e_idx(randi(top_n))).UAV_pos, N_UAV, 2);
        else
            leader_G1 = bestUAV;
            leader_G2 = bestUAV;
            leader_G3 = bestUAV;
        end

        for i = 1:params.K
            cand_i = squeeze(candidates(i, :, :));
            if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                cand_i = cand_i';
            end

            for uav_idx = 1:N_UAV
                if g == 1
                    subpop_best_uav = leader_G1(uav_idx, :);
                elseif g == 2
                    subpop_best_uav = leader_G2(uav_idx, :);
                else
                    subpop_best_uav = leader_G3(uav_idx, :);
                end

                cap_eff = capturability_g(g) * (0.65 + 0.35 * (1 - phi_t));  % 保障Pareto leader引导力，维持前沿覆盖
                pos = goaTurn(cand_i(uav_idx, :), subpop_best_uav, cap_eff, t);
                pos = projectToFeasiblePosition(pos, cand_i(uav_idx, :), cand_i, uav_idx, RRH, N_RRH, N_UAV, params, Ub, Lb);
                candidates(i, uav_idx, :) = pos(:)';
            end

            candidate_pos = squeeze(candidates(i, :, :));
            if size(candidate_pos, 1) == 2 && size(candidate_pos, 2) == N_UAV
                candidate_pos = candidate_pos';
            end
            if ~isFeasibleCandidate(candidate_pos, RRH, N_RRH, N_UAV, params)
                cand_init_i = squeeze(candidates_init(i, :, :));
                if size(cand_init_i, 1) == 2 && size(cand_init_i, 2) == N_UAV
                    cand_init_i = cand_init_i';
                end
                candidate_pos = cand_init_i;
            end
            candidates(i, :, :) = reshape(candidate_pos, 1, N_UAV, 2);
        end

        mem_matrix{g} = updateMemory(mem_matrix{g}, candidates, User, priorities, ...
            E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
    end

    local_mus = zeros(3, N_UAV, 2);
    for g = 1:3
        for uav_idx = 1:N_UAV
            local_mus(g, uav_idx, :) = mean(squeeze(mem_matrix{g}(:, uav_idx, :)), 1);
        end
        subpops{g} = updateSubpopPV(subpops{g}, mem_matrix{g}, squeeze(local_mus(g, :, :)), params.subpop_params, g, iter, params.FES_max, N_UAV);
    end

    [curr_fit_best, curr_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
        User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

    curr_curve(iter) = curr_fit_best;

    if curr_fit_best > scalar_best_fit + 1e-6
        scalar_best_fit = curr_fit_best;
        bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
            E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
        [~, weighted_best, ~, ~] = calcFitness(bestUAV, User, priorities, ...
            E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
    else
        weighted_best = weighted_best_curve(iter - 1);
    end

    [scalar_util, scalar_lat, scalar_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
    cg_curve(iter) = scalar_util;
    weighted_best_curve(iter) = weighted_best;
    energy_consumption(iter) = curr_energy;

    curr_subpop_fits = zeros(3, 1);
    for g = 1:3
        subpop_fits = zeros(1, size(mem_matrix{g}, 1));
        for i = 1:size(mem_matrix{g}, 1)
            candidate = squeeze(mem_matrix{g}(i, :, :));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end
        curr_subpop_fits(g) = max(subpop_fits);
    end

    improvement_threshold_subpop = 1e-6;
    for g = 1:3
        if curr_subpop_fits(g) > prev_fits(g) + improvement_threshold_subpop
            stagnation_counter(g) = 0;
            prev_fits(g) = curr_subpop_fits(g);
        else
            stagnation_counter(g) = stagnation_counter(g) + 1;
        end

        if iter > 20 && stagnation_counter(g) >= 20
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

    E_remaining = params.E_max * ones(N_UAV, 1);
    E_remaining_history(iter, :) = E_remaining';

    for g = 1:3
        subpops{g}.prev_mu = subpops{g}.mu;
    end

    if mod(iter, 50) == 0 || iter == params.FES_max
        fprintf('迭代 %d/%d, 综合得分: %.4f | 当前标量最优解真实效用: %.1f | 时延: %.2f s | 能耗: %.2f J | 档案解数量: %d\n', ...
            iter, params.FES_max, scalar_best_fit, scalar_util, scalar_lat, scalar_nrg, length(pareto_archive));
    end

    pareto_updated_this_iter = false;

    do_pv_mix = params.enable_pv_interpolation && params.enable_multi_subpop && ...
        iter >= 20 && mod(iter, params.pv_interpolation_interval) == 0 && ...
        length(pareto_archive) >= params.pv_interpolation_min_archive && ...
        rand < pv_accept;
    if do_pv_mix
        mixed_candidates = pvInterpolationExchange(subpops, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU, User, priorities, params);
        for mc = 1:length(mixed_candidates)
            cand_pos = mixed_candidates(mc).UAV_pos;
            [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(cand_pos, User, priorities, params);
            [pareto_archive, is_updated] = updateParetoArchive3D(pareto_archive, cand_pos, cand_util, cand_lat, cand_nrg);
            if is_updated
                pareto_updated_this_iter = true;
            end
        end
    end

    for g = 1:3
        for i = 1:size(mem_matrix{g}, 1)
            candidate = squeeze(mem_matrix{g}(i, :, :));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(candidate, User, priorities, params);
            [pareto_archive, is_updated] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
            if is_updated
                pareto_updated_this_iter = true;
            end
        end
    end

    if params.enable_smart_stop && params.enable_early_stop && iter > 100
        if pareto_updated_this_iter
            mo_stagnation_counter = 0;
        else
            mo_stagnation_counter = mo_stagnation_counter + 1;
        end

        if mo_stagnation_counter >= 25
            fprintf('\n[MO-Aware Smart Stop] 迭代 %d/%d: 连续 %d 代未能扩展 Pareto 前沿\n', ...
                iter, params.FES_max, 25);
            actual_iter = iter;
            break;
        end
    end

    actual_iter = iter;
end

center_point = repmat([500, 500], N_UAV, 1);
fly_distances = sqrt(sum((bestUAV - center_point).^2, 2));
fly_energy_costs = params.k_move * fly_distances;
final_E_remaining = max(0, params.E_max - fly_energy_costs');

[final_util, final_lat, final_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
[pareto_archive, ~] = updateParetoArchive3D(pareto_archive, bestUAV, final_util, final_lat, final_nrg);

best_scalar_solution = buildSolutionStruct(bestUAV, scalar_best_fit, User, priorities, params, 'scalar_best');
best_utility_solution = selectMaxUtilitySolution(pareto_archive, bestUAV, scalar_best_fit, User, priorities, params);
best_fit = best_utility_solution.Utility;

fprintf('\n========== 算法完成 ==========\n');
if actual_iter < params.FES_max
    fprintf('提前停止于迭代 %d/%d\n', actual_iter, params.FES_max);
else
    fprintf('达到最大迭代次数 %d\n', params.FES_max);
end
fprintf('最终结果：\n');
fprintf('  1. 历史最佳标量得分：%.4f\n', best_scalar_solution.ScalarFitness);
fprintf('  2. 当前标量最优解真实效用：%.1f | 时延: %.4f s | 能耗: %.2f J\n', ...
    best_scalar_solution.Utility, best_scalar_solution.Latency, best_scalar_solution.Energy);
fprintf('  3. 档案最大 Utility 解：%.1f | 时延: %.4f s | 能耗: %.2f J\n', ...
    best_utility_solution.Utility, best_utility_solution.Latency, best_utility_solution.Energy);
fprintf('  4. Pareto 存档非支配解数量：%d 个\n', length(pareto_archive));
fprintf('============================\n');

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
