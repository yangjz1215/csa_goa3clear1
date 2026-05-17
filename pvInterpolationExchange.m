function mixed_candidates = pvInterpolationExchange(subpops, N_UAV, Ub, Lb, RRH, D_UU, D_RU, User, priorities, params)
%PVINTERPOLATIONEXCHANGE 子群PV凸组合采样：每对子群产生贪心候选
%   贪心候选在α网格上搜索最优标量J
%   接口不变：返回struct数组，含UAV_pos/alpha/pair/source字段

    n_subpops = length(subpops);
    if n_subpops < 2
        mixed_candidates = [];
        return;
    end

    pairs = [2, 3; 1, 3; 1, 2];
    mixed_candidates = [];

    N_User = size(User, 1);
    maxL = N_User * params.max_latency + 1e-9;
    maxE = N_UAV * params.E_max + 1e-9;
    maxU = sum(priorities) + 1e-9;
    if isfield(params, 'G_weights') && numel(params.G_weights) >= 3
        wu = params.G_weights(1);
        wl = params.G_weights(2);
        we = params.G_weights(3);
    else
        wu = 0.4;
        wl = 0.3;
        we = 0.3;
    end

    alpha_grid = [0, 0.5, 1];

    for p = 1:size(pairs, 1)
        idx_a = pairs(p, 1);
        idx_b = pairs(p, 2);

        mu_a = subpops{idx_a}.mu;
        sigma_a = subpops{idx_a}.sigma;
        mu_b = subpops{idx_b}.mu;
        sigma_b = subpops{idx_b}.sigma;

        % --- 候选A：α网格贪心搜索 ---
        best_j = -inf;
        best_pos = [];
        best_alpha = 0.5;

        for ka = 1:numel(alpha_grid)
            alpha = alpha_grid(ka);
            [uav_pos, valid] = mixAndSample(mu_a, sigma_a, mu_b, sigma_b, alpha, N_UAV, Ub, Lb, D_UU);
            if ~valid, continue; end

            [util, lat, nrg] = calcMEC_Objectives(uav_pos, User, priorities, params);
            J = wu * (util / maxU) - wl * (lat / maxL) - we * (nrg / maxE);

            if J > best_j
                best_j = J;
                best_pos = uav_pos;
                best_alpha = alpha;
            end
        end

        if ~isempty(best_pos)
            cand = struct();
            cand.UAV_pos = best_pos;
            cand.alpha = best_alpha;
            cand.pair = [idx_a, idx_b];
            cand.source = 'pv_greedy';
            mixed_candidates = [mixed_candidates; cand];
        end

    end
end

function [uav_pos, valid] = mixAndSample(mu_a, sigma_a, mu_b, sigma_b, alpha, N_UAV, Ub, Lb, D_UU)
%MIXANDSAMPLE 对两个子群PV做凸混合后采样一个UAV部署方案
    mu_mix = alpha * mu_a + (1 - alpha) * mu_b;
    sigma_mix = alpha * sigma_a + (1 - alpha) * sigma_b;

    uav_pos = zeros(N_UAV, 2);
    for uav_idx = 1:N_UAV
        pos = mu_mix(uav_idx, :) + sigma_mix(uav_idx, :) .* randn(1, 2);
        uav_pos(uav_idx, :) = max(Lb, min(Ub, pos));
    end
    uav_pos = enforceUAVDistanceLocal(uav_pos, D_UU);
    valid = true;
end

function pos = enforceUAVDistanceLocal(pos, D_UU)
    N = size(pos, 1);
    max_iter = 10;  % 收敛通常2-3轮，上限10足够
    for iter = 1:max_iter
        moved = false;
        for i = 1:N
            for j = i+1:N
                dist = norm(pos(i, :) - pos(j, :));
                if dist < D_UU && dist > 0
                    dir = randn(1, 2);
                    if norm(dir) < 1e-6
                        dir = [1, 0];
                    end
                    dir = dir / norm(dir);
                    delta = (D_UU * 1.1 - dist) / 2;
                    pos(i, :) = pos(i, :) - dir * delta;
                    pos(j, :) = pos(j, :) + dir * delta;
                    moved = true;
                end
            end
        end
        if ~moved
            break;
        end
    end
    % 兜底：确保严格满足最小距离
    for i = 1:N
        for j = i+1:N
            dist = norm(pos(i, :) - pos(j, :));
            if dist < D_UU && dist > 0
                dir = (pos(i, :) - pos(j, :)) / (dist + 1e-6);
                if norm(dir) < 1e-6
                    dir = [1, 0];
                end
                dir = dir / norm(dir);
                mid_point = (pos(i, :) + pos(j, :)) / 2;
                pos(i, :) = mid_point + dir * (D_UU * 1.1 / 2);
                pos(j, :) = mid_point - dir * (D_UU * 1.1 / 2);
            end
        end
    end
end
