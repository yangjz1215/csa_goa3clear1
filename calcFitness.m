function [fit, obj_utility, obj_latency, obj_energy] = calcFitness(uav_pos, User, priorities, E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params)
    uav_pos = reshape(uav_pos, N_UAV, 2);
    N_User = size(User, 1);

    obj_utility = 0;
    obj_latency = 0;
    obj_energy = 0;

    [utility, latency, energy] = calcMEC_Objectives(uav_pos, User, priorities, params);

    obj_utility = utility;
    obj_latency = latency;
    obj_energy = energy;

    max_utility = sum(priorities);
    max_latency = N_User * params.max_latency;
    if isfield(params, 'energy_norm_max') && params.energy_norm_max > 0
        max_energy = params.energy_norm_max;
    else
        max_energy = E_max * N_UAV;
    end

    w1_def = [0.70, 0.15, 0.15];
    w2_def = [0.30, 0.50, 0.20];
    w3_def = [0.20, 0.15, 0.65];

    if isfield(params, 'test_weights')
        w1 = params.test_weights(g, 1);
        w2 = params.test_weights(g, 2);
        w3 = params.test_weights(g, 3);
    else
        w1 = w1_def(g); w2 = w2_def(g); w3 = w3_def(g);
    end

    norm_util = utility / max_utility;
    norm_lat = latency / max_latency;
    norm_energy = energy / max_energy;
    
    % 基础的多目标适应度 (最大化逻辑：效用越大越好，延迟/能耗越小越好)
    base_fit = w1 * norm_util + w2 * (1.0 - norm_lat) + w3 * (1.0 - norm_energy);
    
    % 防扎堆惩罚与分散探索奖励
    uav_dispersion = 0;
    distance_penalty = 0;
    
    if N_UAV > 1 && isfield(params, 'D_UU')
        min_distances = zeros(N_UAV, 1);
        for i = 1:N_UAV
            dists = sqrt(sum((uav_pos - repmat(uav_pos(i,:), N_UAV, 1)).^2, 2));
            dists(i) = Inf;
            min_distances(i) = min(dists);
            
            for j = i+1:N_UAV
                dist_ij = norm(uav_pos(i,:) - uav_pos(j,:));
                if dist_ij < params.D_UU
                    distance_penalty = distance_penalty + 5.0 * ((params.D_UU - dist_ij) / params.D_UU)^2;
                end
            end
        end
        uav_dispersion = mean(min_distances) / (params.D_UU * 2);
        uav_dispersion = min(1, uav_dispersion);
    end
    
    % 最终适应度 = 基础适应度 + 分散度奖励 - 扎堆惩罚
    fit = base_fit + (0.05 * uav_dispersion) - distance_penalty;
    fit = max(0, double(fit(1)));
end
