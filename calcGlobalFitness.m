function [global_fit, total_energy] = calcGlobalFitness(mem_matrix, G_weights, User, priorities, E_remaining, E_max, k_move, subpop_params, N_UAV, cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params)
    % 使用传入的G_weights，而不是硬编码的权重
    weights = G_weights;
    n_subpops = length(G_weights);  % 动态获取子种群数量
    fits = zeros(1, n_subpops);
    success_best = zeros(1, n_subpops);
    
    for g = 1:n_subpops
        if nargin >= 12 && ~isempty(capturability_g) && length(capturability_g) >= g
            cap = capturability_g(g);
        else
            cap = 0.5;
        end
        
        subpop_size = size(mem_matrix{g},1);
        subpop_fits = zeros(1, subpop_size);
        subpop_success = zeros(1, subpop_size);
        for i = 1:size(mem_matrix{g},1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, cap, N_RRH, RRH_type, UAV_type, params);
        end
        fits(g) = max(subpop_fits);
        success_best(g) = max(subpop_success);
    end
    
    if sum(success_best) > 0
        success_weights = success_best / sum(success_best);
        weights = 0.5 * weights + 0.5 * success_weights;
        weights = weights / sum(weights);
    end
    
    global_fit = sum(fits .* weights);
    global_fit = max(0, global_fit);

    % 计算当前最优解的飞行能耗（基于位置）
    [~, best_g] = max(fits);
    [~, best_local_idx] = max(subpop_fits);

    best_candidate = mem_matrix{best_g}(best_local_idx,:,:);
    best_candidate = reshape(best_candidate, N_UAV, 2);

    center_point = repmat([500, 500], N_UAV, 1);
    fly_distances = sqrt(sum((best_candidate - center_point).^2, 2));
    total_energy = sum(k_move * fly_distances);
end