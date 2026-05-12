function mixed_candidates = pvInterpolationExchange(subpops, N_UAV, Ub, Lb, RRH, D_UU, D_RU, n_alphas)
    if nargin < 9
        n_alphas = 3;
    end

    n_subpops = length(subpops);
    if n_subpops < 2
        mixed_candidates = [];
        return;
    end

    % Prefer G2–G3 (时延–能耗) and G1–G3 (效用–能耗) before G1–G2，利于折中区域与 IGD
    pairs = [2, 3; 1, 3; 1, 2];
    mixed_candidates = [];

    for p = 1:size(pairs, 1)
        idx_a = pairs(p, 1);
        idx_b = pairs(p, 2);

        mu_a = subpops{idx_a}.mu;
        sigma_a = subpops{idx_a}.sigma;
        mu_b = subpops{idx_b}.mu;
        sigma_b = subpops{idx_b}.sigma;

        for k = 1:n_alphas
            alpha = rand();

            mu_mix = alpha * mu_a + (1 - alpha) * mu_b;
            sigma_mix = alpha * sigma_a + (1 - alpha) * sigma_b;

            uav_pos = zeros(N_UAV, 2);
            for uav_idx = 1:N_UAV
                pos = mu_mix(uav_idx, :) + sigma_mix(uav_idx, :) .* randn(1, 2);
                uav_pos(uav_idx, :) = max(Lb, min(Ub, pos));
            end

            uav_pos = enforceUAVDistanceLocal(uav_pos, D_UU);

            candidate = struct();
            candidate.UAV_pos = uav_pos;
            candidate.alpha = alpha;
            candidate.pair = [idx_a, idx_b];
            candidate.source = 'pv_interpolation';

            mixed_candidates = [mixed_candidates; candidate];
        end
    end
end

function pos = enforceUAVDistanceLocal(pos, D_UU)
    N = size(pos, 1);
    max_iter = 20;
    for iter = 1:max_iter
        moved = false;
        for i = 1:N
            for j = i+1:N
                dist = norm(pos(i,:) - pos(j,:));
                if dist < D_UU && dist > 0
                    dir = randn(1, 2);
                    if norm(dir) < 1e-6
                        dir = [1, 0];
                    end
                    dir = dir / norm(dir);
                    delta = (D_UU * 1.1 - dist) / 2;
                    pos(i,:) = pos(i,:) - dir * delta;
                    pos(j,:) = pos(j,:) + dir * delta;
                    moved = true;
                end
            end
        end
        if ~moved
            break;
        end
    end
    for i = 1:N
        for j = i+1:N
            dist = norm(pos(i,:) - pos(j,:));
            if dist < D_UU && dist > 0
                dir = (pos(i,:) - pos(j,:)) / (dist + 1e-6);
                if norm(dir) < 1e-6
                    dir = [1, 0];
                end
                dir = dir / norm(dir);
                mid_point = (pos(i,:) + pos(j,:)) / 2;
                pos(i,:) = mid_point + dir * (D_UU * 1.1 / 2);
                pos(j,:) = mid_point - dir * (D_UU * 1.1 / 2);
            end
        end
    end
end