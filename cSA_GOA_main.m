function [best_fit, bestUAV, cg_curve, energy_consumption, E_remaining_history, final_E_remaining, curr_curve, actual_iter, weighted_best_curve, pareto_archive, best_scalar_solution, best_utility_solution] = ...
    cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)

if ~isfield(params, 'enable_early_stop')
    params.enable_early_stop = true;
end
if ~isfield(params, 'enable_smart_stop')
    params.enable_smart_stop = true;
end
if ~isfield(params, 'enable_bilevel')
    params.enable_bilevel = true;
end
if ~isfield(params, 'enable_multi_subpop')
    params.enable_multi_subpop = true;
end
if ~isfield(params, 'mem_quota_m') || isempty(params.mem_quota_m)
    params.mem_quota_m = 2;
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

fprintf('初始化完成：综合适应度=%.4f | 真实效用(优先级和)=%.1f | 时延=%.2fs | 能耗=%.1f J\n', ...
    scalar_best_fit, init_util, init_lat, init_nrg);

mo_stagnation_counter = 0;
actual_iter = params.FES_max;

% 自适应权重初始化
adaptive_weights = [0.70, 0.15, 0.15;
                    0.30, 0.50, 0.20;
                    0.05, 0.10, 0.85];
if ~isfield(params, 'enable_adaptive_weight')
    params.enable_adaptive_weight = true;
end
params.test_weights = adaptive_weights;

% 适应度缓存：避免 updateMemory 和 calcGlobalFitness 重复评估
mem_fits_cache = cell(3, 1);
mem_utils_cache = cell(3, 1);
mem_lats_cache = cell(3, 1);
mem_nrgs_cache = cell(3, 1);

for iter = 2:params.FES_max
    t = 1 - iter / params.FES_max;

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

                q_eff = params.subpop_params.q(g);
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

            top_pct = max(0.01, 0.20 - 0.19 * (iter / params.FES_max));
            n_archive = length(pareto_archive);
            top_n = max(3, round(n_archive * top_pct));

            % 拥挤度感知的Leader选择：从top-N中选择目标空间中最稀疏的解
            % 计算top-N中每个解到最近邻的目标空间距离，距离越大=越稀疏=优先当leader
            % 关键修复：动态归一化消除量纲影响（Utility~百级, Latency~十级, Energy~万级）
            [~, sort_u_idx] = sort(arch_utils, 'descend');
            top_u_idx = sort_u_idx(1:top_n);
            top_u_objs = zeros(top_n, 3);
            for j = 1:top_n
                top_u_objs(j, :) = [-arch_utils(top_u_idx(j)), arch_lats(top_u_idx(j)), arch_nrgs(top_u_idx(j))];
            end
            min_objs_u = min(top_u_objs, [], 1); max_objs_u = max(top_u_objs, [], 1);
            range_objs_u = max_objs_u - min_objs_u; range_objs_u(range_objs_u < 1e-9) = 1;
            norm_u_objs = (top_u_objs - repmat(min_objs_u, top_n, 1)) ./ repmat(range_objs_u, top_n, 1);
            nn_dist_u = inf(top_n, 1);
            for j = 1:top_n
                d = sqrt(sum((norm_u_objs - repmat(norm_u_objs(j,:), top_n, 1)).^2, 2));
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
            min_objs_l = min(top_l_objs, [], 1); max_objs_l = max(top_l_objs, [], 1);
            range_objs_l = max_objs_l - min_objs_l; range_objs_l(range_objs_l < 1e-9) = 1;
            norm_l_objs = (top_l_objs - repmat(min_objs_l, top_n, 1)) ./ repmat(range_objs_l, top_n, 1);
            nn_dist_l = inf(top_n, 1);
            for j = 1:top_n
                d = sqrt(sum((norm_l_objs - repmat(norm_l_objs(j,:), top_n, 1)).^2, 2));
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
            min_objs_e = min(top_e_objs, [], 1); max_objs_e = max(top_e_objs, [], 1);
            range_objs_e = max_objs_e - min_objs_e; range_objs_e(range_objs_e < 1e-9) = 1;
            norm_e_objs = (top_e_objs - repmat(min_objs_e, top_n, 1)) ./ repmat(range_objs_e, top_n, 1);
            nn_dist_e = inf(top_n, 1);
            for j = 1:top_n
                d = sqrt(sum((norm_e_objs - repmat(norm_e_objs(j,:), top_n, 1)).^2, 2));
                d(j) = inf;
                nn_dist_e(j) = min(d);
            end
            [~, sp_e] = max(nn_dist_e);
            leader_G3 = reshape(pareto_archive(top_e_idx(sp_e)).UAV_pos, N_UAV, 2);
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

                cap_eff = capturability_g(g);
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

        [mem_matrix{g}, mem_fits_cache{g}, mem_utils_cache{g}, mem_lats_cache{g}, mem_nrgs_cache{g}] = ...
            updateMemory(mem_matrix{g}, candidates, User, priorities, ...
            E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params, ...
            mem_fits_cache{g}, mem_utils_cache{g}, mem_lats_cache{g}, mem_nrgs_cache{g});
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
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params, mem_fits_cache);

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

    % Pareto存档更新：使用updateMemory缓存的目标值，避免重复调用calcMEC_Objectives
    % 性能优化：每5代更新一次Pareto存档（原来每代更新导致运行时间暴增4-5倍）
    % 原因：Pareto更新本身是O(n^2)操作，每代更新3子群x10候选 = 30次插入，触发频繁的支配比较和截断
    if mod(iter, 5) == 0 || iter == params.FES_max
        for g = 1:3
            for i = 1:size(mem_matrix{g}, 1)
                candidate = squeeze(mem_matrix{g}(i, :, :));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                cand_util = mem_utils_cache{g}(i);
                cand_lat = mem_lats_cache{g}(i);
                cand_nrg = mem_nrgs_cache{g}(i);
                [pareto_archive, is_updated] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
                if is_updated
                    pareto_updated_this_iter = true;
                end
            end
        end
    end

    % 自适应权重旋转：每20代根据Pareto存档分布微调子种群权重方向
    % 核心思想：找到前沿上最稀疏的区域，将最近的子种群权重朝该方向旋转
    if params.enable_adaptive_weight && mod(iter, 20) == 0 && iter > 20 && iter < 0.8 * params.FES_max && length(pareto_archive) >= 10
        arch_utils_aw = [pareto_archive.Utility];
        arch_lats_aw = [pareto_archive.Latency];
        arch_nrgs_aw = [pareto_archive.Energy];
        n_arch_aw = length(pareto_archive);

        % 归一化目标值到[0,1]
        range_u_aw = max(arch_utils_aw) - min(arch_utils_aw); if range_u_aw < 1e-9, range_u_aw = 1; end
        range_l_aw = max(arch_lats_aw) - min(arch_lats_aw); if range_l_aw < 1e-9, range_l_aw = 1; end
        range_e_aw = max(arch_nrgs_aw) - min(arch_nrgs_aw); if range_e_aw < 1e-9, range_e_aw = 1; end

        nu_aw = (arch_utils_aw - min(arch_utils_aw)) / range_u_aw;
        nl_aw = (arch_lats_aw - min(arch_lats_aw)) / range_l_aw;
        ne_aw = (arch_nrgs_aw - min(arch_nrgs_aw)) / range_e_aw;

        % 计算每个解到最近邻的目标空间距离（拥挤度指标）
        objs_aw = [nu_aw(:), nl_aw(:), ne_aw(:)];
        crowd_aw = zeros(n_arch_aw, 1);
        for ci = 1:n_arch_aw
            d_aw = sqrt(sum((objs_aw - repmat(objs_aw(ci,:), n_arch_aw, 1)).^2, 2));
            d_aw(ci) = inf;
            crowd_aw(ci) = min(d_aw);
        end

        % 修复：不再朝最稀疏点旋转（极端点≠有价值方向）
        % 改为：计算当前3个子种群权重在目标空间中的覆盖盲区，朝盲区方向旋转
        % 3个理想方向：G1=[1,0,0](Utility), G2=[0,1,0](anti-Latency), G3=[0,0,1](anti-Energy)
        ideal_dirs = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        % 计算每个理想方向被当前子种群覆盖的程度（余弦相似度最大值）
        coverage_score = zeros(3, 1);
        for di = 1:3
            max_cos = 0;
            for gi = 1:3
                cos_val = dot(adaptive_weights(gi,:), ideal_dirs(di,:)) / (norm(adaptive_weights(gi,:)) * norm(ideal_dirs(di,:)) + 1e-9);
                max_cos = max(max_cos, cos_val);
            end
            coverage_score(di) = max_cos;
        end
        % 找到覆盖最弱的理想方向
        [~, weakest_dir] = min(coverage_score);
        % 找到离该理想方向最远的子种群（最适合被旋转过去）
        max_dist = 0; closest_g_aw = 1;
        for gi = 1:3
            dist_val = 1 - dot(adaptive_weights(gi,:), ideal_dirs(weakest_dir,:)) / (norm(adaptive_weights(gi,:)) * norm(ideal_dirs(weakest_dir,:)) + 1e-9);
            if dist_val > max_dist
                max_dist = dist_val;
                closest_g_aw = gi;
            end
        end
        target_dir = ideal_dirs(weakest_dir,:);

        % 旋转强度随迭代递减
        progress_aw = iter / params.FES_max;
        alpha_aw = 0.03 * (1 - progress_aw)^1.5;

        % 将最近子种群的权重朝目标方向旋转
        adaptive_weights(closest_g_aw, :) = (1 - alpha_aw) * adaptive_weights(closest_g_aw, :) + alpha_aw * target_dir;

        % 约束：每个权重在[0.05, 0.90]范围内，权重和为1
        w_aw = adaptive_weights(closest_g_aw, :);
        w_aw = max(0.05, min(0.90, w_aw));
        adaptive_weights(closest_g_aw, :) = w_aw / sum(w_aw);

        % 关键约束：保护子种群专业化——每个子种群的主目标权重不能低于阈值
        % G1主Utility(w1>=0.40), G2主Latency(w2>=0.30), G3主Energy(w3>=0.50)
        min_primary = [0.40, 0.30, 0.50];
        if adaptive_weights(closest_g_aw, closest_g_aw) < min_primary(closest_g_aw)
            adaptive_weights(closest_g_aw, closest_g_aw) = min_primary(closest_g_aw);
            adaptive_weights(closest_g_aw, :) = adaptive_weights(closest_g_aw, :) / sum(adaptive_weights(closest_g_aw, :));
        end

        params.test_weights = adaptive_weights;

        % 关键修复：test_weights驱动calcFitness内的3维目标权重
        % G_weights保持为子种群标量权重(1x3)，用于calcGlobalFitness加权求和，不应被3x3矩阵覆盖
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
