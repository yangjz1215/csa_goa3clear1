function global_best = calcGlobalBest(mem_matrix, G_weights, N_UAV, User, priorities, ...
    E_remaining, E_max, k_move, subpop_params, cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params)
    n_subpops = numel(mem_matrix);

    weights = reshape(G_weights, 1, []);
    if isempty(weights)
        weights = ones(1, n_subpops) / n_subpops;
    elseif numel(weights) ~= n_subpops
        weights = weights(1:min(numel(weights), n_subpops));
        if numel(weights) < n_subpops
            weights(end+1:n_subpops) = 1;
        end
        weights = weights / sum(weights);
    end

    best_pos = zeros(n_subpops, N_UAV, 2);
    best_fits = zeros(n_subpops, 1);
    best_dispersions = zeros(n_subpops, 1);
    best_success = zeros(n_subpops, 1);

    for g = 1:n_subpops
        if nargin >= 12 && ~isempty(capturability_g) && length(capturability_g) >= g
            cap = capturability_g(g);
        else
            cap = 0.5;
        end

        subpop_size = size(mem_matrix{g}, 1);
        subpop_fits = zeros(1, subpop_size);
        subpop_success = zeros(1, subpop_size);
        subpop_dispersions = zeros(1, subpop_size);

        for i = 1:subpop_size
            candidate = squeeze(mem_matrix{g}(i, :, :));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end

            [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, cap, N_RRH, RRH_type, UAV_type, params);
            [~, ~, ~, subpop_success(i)] = calcMEC_Objectives(candidate, User, priorities, params);

            if N_UAV > 1
                min_distances = zeros(N_UAV, 1);
                for uav_idx = 1:N_UAV
                    dists = sqrt(sum((candidate - repmat(candidate(uav_idx, :), N_UAV, 1)).^2, 2));
                    dists(uav_idx) = Inf;
                    min_distances(uav_idx) = min(dists);
                end
                subpop_dispersions(i) = mean(min_distances);
            else
                subpop_dispersions(i) = 0;
            end
        end

        fit_norm = subpop_fits / max(max(subpop_fits), 1e-6);
        success_norm = subpop_success / max(max(subpop_success), 1e-6);
        dispersion_norm = subpop_dispersions / max(max(subpop_dispersions), 1e-6);
        combined_scores = 0.6 * fit_norm + 0.4 * success_norm + 0.3 * dispersion_norm;

        [~, best_idx] = max(combined_scores);
        best_pos(g, :, :) = mem_matrix{g}(best_idx, :, :);
        best_fits(g) = subpop_fits(best_idx);
        best_dispersions(g) = subpop_dispersions(best_idx);
        best_success(g) = subpop_success(best_idx);
    end

    [~, best_disp_idx] = max(best_dispersions);

    if best_fits(best_disp_idx) >= 0.8 * max(best_fits)
        global_best = squeeze(best_pos(best_disp_idx, :, :));
        if size(global_best, 1) == 1 && size(global_best, 2) == N_UAV * 2
            global_best = reshape(global_best, N_UAV, 2);
        end
    else
        adjusted_weights = weights;
        adjusted_weights(best_disp_idx) = adjusted_weights(best_disp_idx) * 1.5;
        adjusted_weights = adjusted_weights / sum(adjusted_weights);

        if sum(best_success) > 0
            success_weights = best_success(:)' / sum(best_success);
            adjusted_weights = 0.5 * adjusted_weights + 0.5 * success_weights;
            adjusted_weights = adjusted_weights / sum(adjusted_weights);
        end

        global_best = zeros(N_UAV, 2);
        for uav_idx = 1:N_UAV
            for g = 1:n_subpops
                subpop_uav_pos = squeeze(best_pos(g, uav_idx, :));
                if size(subpop_uav_pos, 1) > 1
                    subpop_uav_pos = subpop_uav_pos(1, :);
                end
                global_best(uav_idx, :) = global_best(uav_idx, :) + adjusted_weights(g) * subpop_uav_pos;
            end
        end
    end

    if size(global_best, 1) ~= N_UAV || size(global_best, 2) ~= 2
        global_best = reshape(global_best, N_UAV, 2);
    end
end
