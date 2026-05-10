function [utility, latency, energy, success_rate] = calcMEC_Objectives(UAV_pos, User, priorities, params)
    N_User = size(User, 1);
    N_UAV = size(UAV_pos, 1);

    utility = 0;
    latency = 0;
    energy = 0;
    success_rate = 0;

    % Prefer the exact Code_main0-style evaluation chain when task and type
    % metadata are available. This keeps the reported "real utility" aligned
    % with the actual connect/offload constraints instead of a simplified
    % nearest-UAV proxy.
    has_exact_inputs = isfield(params, 'D') && isfield(params, 'C') && isfield(params, 'DT') ...
        && isfield(params, 'RRH') && ~isempty(params.RRH);

    if has_exact_inputs
        [utility, latency, energy, success_rate] = evaluateExactMEC(UAV_pos, User, priorities, params);
    else
        [utility, latency, energy, success_rate] = evaluateSimplifiedMEC(UAV_pos, User, priorities, params);
    end
end

function [utility, latency, energy, success_rate] = evaluateExactMEC(UAV_pos, User, priorities, params)
    N_User = size(User, 1);
    N_UAV = size(UAV_pos, 1);
    RRH = params.RRH;
    N_RRH = size(RRH, 1);

    if isfield(params, 'RRH_type') && numel(params.RRH_type) == N_RRH
        RRH_type = params.RRH_type(:)';
    else
        RRH_type = zeros(1, N_RRH);
    end

    if isfield(params, 'UAV_type') && numel(params.UAV_type) == N_UAV
        UAV_type = params.UAV_type(:)';
    else
        UAV_type = zeros(1, N_UAV);
    end

    H = 50;
    omega0 = 1e6;
    alpha0 = 1.42e-4;
    sigma2 = 3.98e-12;
    f_eRRH = 2e9;
    f_eUAV = 2e9;
    f_BBU = 4e9;
    Ptx = 1;
    PtxR = 10;
    PtxU = 10;
    PtxEU = 5;
    Pho = 100;
    ki = 1e-27;
    k_move = 15;

    if isfield(params, 'alpha0'); alpha0 = params.alpha0; end
    if isfield(params, 'sigma2'); sigma2 = params.sigma2; end
    if isfield(params, 'f_eRRH'); f_eRRH = params.f_eRRH; end
    if isfield(params, 'f_eUAV'); f_eUAV = params.f_eUAV; end
    if isfield(params, 'f_BBU'); f_BBU = params.f_BBU; end
    if isfield(params, 'Ptx'); Ptx = params.Ptx; end
    if isfield(params, 'PtxR'); PtxR = params.PtxR; end
    if isfield(params, 'PtxU'); PtxU = params.PtxU; end
    if isfield(params, 'PtxEU'); PtxEU = params.PtxEU; end
    if isfield(params, 'Pho'); Pho = params.Pho; end
    if isfield(params, 'ki'); ki = params.ki; end
    if isfield(params, 'k_move'); k_move = params.k_move; end

    if isscalar(params.E_max)
        E_max = params.E_max * ones(N_UAV, 1);
    else
        E_max = params.E_max(:);
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_distances = sqrt(sum((UAV_pos - center_point) .^ 2, 2));
    flight_energy = k_move * fly_distances;
    actual_E_remaining = max(0, E_max - flight_energy);

    Connect = UserConnect_EnergyAware(User, N_User, RRH, N_RRH, RRH_type, UAV_pos, N_UAV, UAV_type, ...
        priorities, actual_E_remaining, E_max, params);
    off_decision = UserOffloading_EnergyAware(User, N_User, N_RRH, RRH, RRH_type, UAV_pos, N_UAV, UAV_type, ...
        params.D, params.C, priorities, Connect, actual_E_remaining, E_max, params);

    success_flags = zeros(1, N_User);
    latency_vec = zeros(1, N_User);
    hover_time = zeros(N_UAV, 1);
    relay_energy = 0;
    compute_energy = 0;

    for i = 1:N_User
        if Connect(i) == 0
            continue;
        end

        index = off_decision(i);
        if index == 0
            continue;
        end

        D_user = getPerUserValue(params.D, i);
        C_user = getPerUserValue(params.C, i);
        DT_user = getPerUserValue(params.DT, i);

        if index <= N_RRH
            h0 = alpha0 / sum((RRH(index, :) - User(i, :)) .^ 2);
            r = omega0 * log2(1 + (Ptx * h0) / sigma2);
            t_tx = D_user / r;
            t_ex = C_user / f_eRRH;
            total_latency = t_tx + t_ex;

            if total_latency <= DT_user
                success_flags(i) = 1;
                latency_vec(i) = total_latency;
            end

        elseif index <= N_RRH + N_UAV
            uav_idx = index - N_RRH;
            h0 = alpha0 / sum((UAV_pos(uav_idx, :) - User(i, :)) .^ 2 + H^2);
            tx_power = Ptx;
            r = omega0 * log2(1 + (tx_power * h0) / sigma2);
            t_tx = D_user / r;
            t_ex = C_user / f_eUAV;
            total_latency = t_tx + t_ex;

            if total_latency <= DT_user
                success_flags(i) = 1;
                latency_vec(i) = total_latency;
                hover_time(uav_idx) = max(hover_time(uav_idx), total_latency);
                compute_energy = compute_energy + ki * (f_eUAV ^ 2) * C_user;
            end

        else
            cindex = Connect(i);
            if cindex <= N_RRH
                h0 = alpha0 / sum((RRH(cindex, :) - User(i, :)) .^ 2);
                r1 = omega0 * log2(1 + (Ptx * h0) / sigma2);
                t_tx1 = D_user / r1;

                h01 = alpha0 / sum((RRH(cindex, :)) .^ 2);
                r2 = omega0 * log2(1 + (PtxR * h01) / sigma2);
                t_tx2 = D_user / r2;

                total_latency = t_tx1 + t_tx2 + C_user / f_BBU;
                if total_latency <= DT_user
                    success_flags(i) = 1;
                    latency_vec(i) = total_latency;
                end
            else
                uav_idx = cindex - N_RRH;
                if UAV_type(uav_idx) == 1
                    relay_power = PtxEU;
                else
                    relay_power = PtxU;
                end

                h0 = alpha0 / sum((UAV_pos(uav_idx, :) - User(i, :)) .^ 2 + H^2);
                r1 = omega0 * log2(1 + (Ptx * h0) / sigma2);
                t_tx1 = D_user / r1;

                h01 = alpha0 / sum((UAV_pos(uav_idx, :)) .^ 2 + H^2);
                r2 = omega0 * log2(1 + (relay_power * h01) / sigma2);
                t_tx2 = D_user / r2;

                total_latency = t_tx1 + t_tx2 + C_user / f_BBU;
                if total_latency <= DT_user
                    success_flags(i) = 1;
                    latency_vec(i) = total_latency;
                    hover_time(uav_idx) = max(hover_time(uav_idx), t_tx1 + t_tx2);
                    relay_energy = relay_energy + relay_power * t_tx2;
                end
            end
        end
    end

    utility = sum(priorities(success_flags == 1));
    latency = sum(latency_vec(success_flags == 1));
    energy = sum(flight_energy) + relay_energy + compute_energy + Pho * sum(hover_time);
    success_rate = mean(success_flags == 1);
end

function [utility, latency, energy, success_rate] = evaluateSimplifiedMEC(UAV_pos, User, priorities, params)
    N_User = size(User, 1);
    N_UAV = size(UAV_pos, 1);

    utility = 0;
    latency = 0;
    comp_energy = 0;
    success_count = 0;

    B_available = ones(N_UAV, 1) * params.B_total;
    F_available = ones(N_UAV, 1) * params.F_total;
    B_relay_available = ones(N_UAV, 1) * params.B_total_relay;

    [sorted_prio, sorted_idx] = sort(priorities, 'descend');

    for i = 1:N_User
        u_idx = sorted_idx(i);
        prio = sorted_prio(i);

        task_D_user = params.D(u_idx);
        task_C_user = params.C(u_idx);

        dists = sqrt(sum((UAV_pos - User(u_idx, :)) .^ 2, 2));
        [min_dist, best_uav] = min(dists);
        min_dist = max(min_dist, 10.0);

        if min_dist <= params.cover_radius
            snr = params.P_tx / (min_dist^2 * params.noise);
            spectral_efficiency = log2(1 + snr);

            max_t_trans = params.max_latency * 0.5;
            max_t_comp = params.max_latency * 0.5;

            req_B = task_D_user / (max_t_trans * spectral_efficiency);
            req_F = task_C_user / max_t_comp;

            if B_available(best_uav) >= req_B && F_available(best_uav) >= req_F
                B_available(best_uav) = B_available(best_uav) - req_B;
                F_available(best_uav) = F_available(best_uav) - req_F;

                utility = utility + prio;
                success_count = success_count + 1;

                t_trans = task_D_user / (req_B * spectral_efficiency);
                t_comp = task_C_user / req_F;
                latency = latency + (t_trans + t_comp);

                comp_energy = comp_energy + params.kappa * req_F^2 * task_C_user;
            else
                if isfield(params, 'RRH') && ~isempty(params.RRH) && isfield(params, 'PtxU') ...
                        && isfield(params, 'B_total_relay') && isfield(params, 'f_BBU')
                    dists_rrh = sqrt(sum((params.RRH - UAV_pos(best_uav, :)) .^ 2, 2));
                    [min_dist_rrh, ~] = min(dists_rrh);
                    min_dist_rrh = max(min_dist_rrh, 10.0);

                    if min_dist_rrh <= 300
                        req_B_access = task_D_user / ((params.max_latency * 0.4) * spectral_efficiency);

                        if B_available(best_uav) >= req_B_access
                            snr_relay = params.PtxU / (min_dist_rrh^2 * params.noise);
                            spectral_eff_relay = log2(1 + snr_relay);

                            req_B_relay = task_D_user / ((params.max_latency * 0.4) * spectral_eff_relay);

                            if B_relay_available(best_uav) >= req_B_relay
                                B_available(best_uav) = B_available(best_uav) - req_B_access;
                                B_relay_available(best_uav) = B_relay_available(best_uav) - req_B_relay;

                                t_trans1 = task_D_user / (req_B_access * spectral_efficiency);
                                t_trans2 = task_D_user / (req_B_relay * spectral_eff_relay);
                                t_comp_bbu = task_C_user / params.f_BBU;
                                total_relay_latency = t_trans1 + t_trans2 + t_comp_bbu;

                                if total_relay_latency <= params.max_latency
                                    utility = utility + prio;
                                    latency = latency + total_relay_latency;
                                    comp_energy = comp_energy + (params.PtxU * t_trans2);
                                    success_count = success_count + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    fly_dist = sqrt(sum((UAV_pos - repmat([500, 500], N_UAV, 1)) .^ 2, 2));
    flight_energy = sum(fly_dist * params.k_move);

    energy = flight_energy + comp_energy;
    success_rate = success_count / N_User;
end

function value = getPerUserValue(field_value, idx)
    if isscalar(field_value)
        value = field_value;
    else
        value = field_value(idx);
    end
end
