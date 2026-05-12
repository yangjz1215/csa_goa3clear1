function phi = computePhasePhi(iter, FES_max, bestUAV, User, priorities, params, RRH)
%COMPUTEPHASEPHI 灾后双层场景序参量 φ∈[0,1]：覆盖缺口 + 内层(时延/能耗)压力 + 迭代进程
%   φ 高：灾情/网络压力大或早期 → 偏探索；φ 低：后期稳定 → 偏折中与精修。
%   仅用现有 MEC 输出与覆盖比，不新增约束、不改 calcMEC_Objectives 内部。

    N_User = size(User, 1);
    N_UAV = size(bestUAV, 1);
    progress = iter / max(FES_max, 1);

    cov = localCoverageRatio(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);
    cov_gap = max(0, min(1, 1 - cov));

    [u, lat, nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
    maxL = N_User * params.max_latency + 1e-9;
    maxE = N_UAV * params.E_max + 1e-9;
    lat_n = min(1, lat / maxL);
    nrg_n = min(1, nrg / maxE);
    inner = 0.5 * lat_n + 0.5 * nrg_n;

    if isfield(params, 'phase_w_progress') && ~isempty(params.phase_w_progress)
        w1 = params.phase_w_progress;
    else
        w1 = 0.28;
    end
    if isfield(params, 'phase_w_cov') && ~isempty(params.phase_w_cov)
        w2 = params.phase_w_cov;
    else
        w2 = 0.42;
    end
    if isfield(params, 'phase_w_inner') && ~isempty(params.phase_w_inner)
        w3 = params.phase_w_inner;
    else
        w3 = 0.30;
    end
    s = w1 + w2 + w3;
    w1 = w1 / s;
    w2 = w2 / s;
    w3 = w3 / s;

    phi = w1 * (1 - progress) + w2 * cov_gap + w3 * inner;
    phi = min(1, max(0, phi));
end

function cov_ratio = localCoverageRatio(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos, 1)
        dists_uav = sqrt(sum((UAV_pos - User_pos(i, :)).^2, 2));
        covered_by_uav = any(dists_uav <= UAV_radius);
        if size(RRH, 1) > 0
            dists_rrh = sqrt(sum((RRH - User_pos(i, :)).^2, 2));
            covered_by_rrh = any(dists_rrh <= RRH_radius);
        else
            covered_by_rrh = false;
        end
        if covered_by_uav || covered_by_rrh
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos, 1);
end
