function subpops = initSubpopulations(N_UAV, User, RRH, priorities, subpop_params, Ub, Lb, cover_radius, D_RU)
    subpops = cell(3,1);
    D_UU = 10;

    start_point = (Ub + Lb) / 2;

    g1_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    subpops{1}.mu = g1_pos;
    subpops{1}.sigma = repmat(subpop_params.sigma0(1) * 0.3, N_UAV, 2);
    subpops{1}.prev_mu = subpops{1}.mu;
    subpops{1}.id = 1;

    g2_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    subpops{2}.mu = g2_pos;
    subpops{2}.sigma = repmat(subpop_params.sigma0(2) * 0.3, N_UAV, 2);
    subpops{2}.prev_mu = subpops{2}.mu;
    subpops{2}.id = 2;

    if size(RRH,1) > 0
        ref_rrh_idx = randi(size(RRH, 1));
        ref_rrh = RRH(ref_rrh_idx, :);
        g3_pos = generateCenteredUAV(N_UAV, ref_rrh, RRH, D_UU, D_RU, Ub, Lb, 2);
    else
        g3_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    end
    subpops{3}.mu = g3_pos;
    subpops{3}.sigma = repmat(subpop_params.sigma0(3) * 0.3, N_UAV, 2);
    subpops{3}.prev_mu = subpops{3}.mu;
    subpops{3}.id = 3;
end

function pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, jitter)
    pos = repmat(start_point, N_UAV, 1);
    max_attempts = 100;
    for i = 1:N_UAV
        valid = false;
        attempts = 0;
        while ~valid && attempts < max_attempts
            candidate = start_point + jitter * randn(1, 2);
            candidate = max(Lb, min(Ub, candidate));
            if i == 1
                valid = true;
            else
                dist_to_others = sqrt(sum((pos(1:i-1,:) - repmat(candidate, i-1, 1)).^2, 2));
                dist_to_rrh = sqrt(sum((RRH - candidate).^2, 2));
                valid = all(dist_to_others >= D_UU) && all(dist_to_rrh >= D_RU);
            end
            attempts = attempts + 1;
        end
        if valid
            pos(i,:) = candidate;
        else
            pos(i,:) = start_point;
        end
    end
end