function [global_fit, total_energy] = calcGlobalFitness(mem_matrix, G_weights, User, priorities, E_remaining, E_max, k_move, subpop_params, N_UAV, cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params, cached_fits_cell)
%CALCGLOBALFITNESS 计算全局加权适应度
%   优化：接受 cached_fits_cell（每个子群的预计算适应度），跳过重复评估

    weights = G_weights;
    n_subpops = length(G_weights);
    fits = zeros(1, n_subpops);
    best_local_indices = zeros(1, n_subpops);

    has_cache = (nargin >= 17 && ~isempty(cached_fits_cell));

    for g = 1:n_subpops
        if nargin >= 12 && ~isempty(capturability_g) && length(capturability_g) >= g
            cap = capturability_g(g);
        else
            cap = 0.5;
        end

        subpop_size = size(mem_matrix{g},1);

        % 如果有缓存且维度匹配，直接使用
        use_cache = has_cache && length(cached_fits_cell) >= g && ~isempty(cached_fits_cell{g}) && length(cached_fits_cell{g}) == subpop_size;

        if use_cache
            subpop_fits = cached_fits_cell{g};
        else
            subpop_fits = zeros(1, subpop_size);
            for i = 1:subpop_size
                candidate = squeeze(mem_matrix{g}(i,:,:));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                    E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, cap, N_RRH, RRH_type, UAV_type, params);
            end
        end
        [fits(g), best_local_indices(g)] = max(subpop_fits);
    end

    global_fit = sum(fits .* weights);
    global_fit = max(0, global_fit);

    [~, best_g] = max(fits);
    best_local_idx = best_local_indices(best_g);

    best_candidate = mem_matrix{best_g}(best_local_idx,:,:);
    best_candidate = reshape(best_candidate, N_UAV, 2);

    center_point = repmat([500, 500], N_UAV, 1);
    fly_distances = sqrt(sum((best_candidate - center_point).^2, 2));
    total_energy = sum(k_move * fly_distances);
end
