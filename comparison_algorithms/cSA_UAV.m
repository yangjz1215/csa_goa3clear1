function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
    end

    if isfield(params, 'K')
        pop_size = params.K * 3;
    else
        pop_size = 120;
    end

    if ~isfield(params, 'enable_bilevel'); params.enable_bilevel = true; end
    if ~isfield(params, 'B_total'); params.B_total = 20e6; end
    if ~isfield(params, 'F_total'); params.F_total = 10e9; end
    if ~isfield(params, 'max_latency'); params.max_latency = 1.0; end
    if ~isfield(params, 'kappa'); params.kappa = 1e-27; end
    if ~isfield(params, 'P_tx'); params.P_tx = 1; end
    if ~isfield(params, 'noise'); params.noise = 1e-13; end
    if ~isfield(params, 'RRH_radius'); params.RRH_radius = 120; end
    params.RRH = RRH;

    n_vars = N_UAV * 2;
    center_point = [500, 500];
    jitter = 10;
    Np = pop_size;

    population = zeros(pop_size, n_vars);
    fitness_pop = zeros(pop_size, 1);
    mu = zeros(pop_size, n_vars);
    sigma = ones(pop_size, n_vars);

    E_remaining = params.E_max * ones(N_UAV, 1);
    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

    for i = 1:pop_size
        for j = 1:N_UAV
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            population(i, (j-1)*2+1) = max(Lb(1), min(Ub(1), init_x));
            population(i, (j-1)*2+2) = max(Lb(2), min(Ub(2), init_y));
        end

        uav_pos = reshape(population(i, :), N_UAV, 2);
        if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
            uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
            population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
        end

        center_point_rep = repmat([500, 500], N_UAV, 1);
        fly_dist = sqrt(sum((uav_pos - center_point_rep).^2, 2));
        fly_energy = params.k_move * fly_dist;
        E_curr = max(0, params.E_max - fly_energy);

        [fitness, ~, ~, ~] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
            params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
        fitness_pop(i) = fitness;

        [cand_util, cand_lat, cand_nrg, ~] = calcMEC_Objectives(uav_pos, User, priorities, params);
        pareto_archive = updateParetoArchive3D(pareto_archive, uav_pos, cand_util, cand_lat, cand_nrg);
    end

    [best_fit, best_idx] = max(fitness_pop);
    bestUAV = reshape(population(best_idx, :), N_UAV, 2);
    center_point_rep = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point_rep).^2, 2));
    best_energy = sum(params.k_move * fly_dist);

    cg_curve = zeros(1, max_iter);
    [init_real_utility, ~, ~, ~] = calcMEC_Objectives(bestUAV, User, priorities, params);
    cg_curve(1) = init_real_utility;

    for iter = 2:max_iter
        a = 2 - iter * (2 / max_iter);

        for i = 1:pop_size
            new_pos = zeros(1, n_vars);
            for j = 1:n_vars
                r1 = 2 * a * rand() - a;
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                new_pos(j) = new_pos(j) + r1 * sin(r2) * (r3 * bestUAV(j) - new_pos(j));
            end

            new_pos = max(Lb(1), min(Ub(1), new_pos));
            uav_pos = reshape(new_pos, N_UAV, 2);
            if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
                uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
                new_pos = reshape(uav_pos, 1, N_UAV * 2);
            end

            winner = new_pos;
            loser = population(i, :);

            center_point_rep = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((reshape(new_pos, N_UAV, 2) - center_point_rep).^2, 2));
            fly_energy = params.k_move * fly_dist;
            E_curr = max(0, params.E_max - fly_energy);

            [candidate_fit, ~, ~, ~] = calcFitness(reshape(new_pos, N_UAV, 2), User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

            if candidate_fit > fitness_pop(i)
                winner = new_pos;
                loser = population(i, :);
                fitness_pop(i) = candidate_fit;
                population(i, :) = new_pos;

                if candidate_fit > best_fit
                    best_fit = candidate_fit;
                    bestUAV = reshape(new_pos, N_UAV, 2);
                    fly_dist = sqrt(sum((bestUAV - center_point_rep).^2, 2));
                    best_energy = sum(params.k_move * fly_dist);
                end
            else
                winner = population(i, :);
                loser = new_pos;
            end

            mu(i, :) = mu(i, :) + (1 / Np) * (winner - loser);
            var_update = sigma(i, :).^2 + mu(i, :).^2 - (mu(i, :) + (1 / Np) * (winner - loser)).^2 ...
                + (1 / Np) * (winner.^2 - loser.^2);
            sigma(i, :) = sqrt(max(var_update, 1e-6));

            [cand_util, cand_lat, cand_nrg, ~] = calcMEC_Objectives(reshape(population(i, :), N_UAV, 2), User, priorities, params);
            pareto_archive = updateParetoArchive3D(pareto_archive, reshape(population(i, :), N_UAV, 2), cand_util, cand_lat, cand_nrg);
        end

        [real_utility, ~, ~, ~] = calcMEC_Objectives(bestUAV, User, priorities, params);
        cg_curve(iter) = real_utility;

        if mod(iter, 50) == 0
            fprintf('cSA iter %d/%d, Best fitness: %.4f\n', iter, max_iter, best_fit);
        end
    end

    if ~isempty(pareto_archive) && length(pareto_archive) >= 1
        arch_utils = [pareto_archive.Utility];
        [~, max_u_idx] = max(arch_utils);
        bestUAV = pareto_archive(max_u_idx).UAV_pos;
        best_fit = pareto_archive(max_u_idx).Utility;
        best_energy = pareto_archive(max_u_idx).Energy;
    else
        center_point_rep = repmat([500, 500], N_UAV, 1);
        fly_dist = sqrt(sum((bestUAV - center_point_rep).^2, 2));
        fly_energy = params.k_move * fly_dist;
        E_curr = max(0, params.E_max - fly_energy);
        [best_fit, best_util, ~, best_nrg] = calcFitness(bestUAV, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
            params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
        best_fit = best_util;
        best_energy = best_nrg;
    end
end

function feasible = checkConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    feasible = true;

    for i = 1:N_UAV
        for j = i+1:N_UAV
            dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
            if dist < D_UU
                feasible = false;
                return;
            end
        end
    end

    for i = 1:N_UAV
        for j = 1:size(RRH, 1)
            dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
            if dist < D_RU
                feasible = false;
                return;
            end
        end
    end
end

function uav_pos = enforceConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    max_iter = 100;

    for iter = 1:max_iter
        changed = false;

        for i = 1:N_UAV
            for j = i+1:N_UAV
                dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
                if dist < D_UU && dist > 0
                    direction = (uav_pos(i, :) - uav_pos(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_UU - dist) * 0.5;
                    uav_pos(j, :) = uav_pos(j, :) - direction * (D_UU - dist) * 0.5;
                    changed = true;
                end
            end
        end

        for i = 1:N_UAV
            for j = 1:size(RRH, 1)
                dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
                if dist < D_RU && dist > 0
                    direction = (uav_pos(i, :) - RRH(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_RU - dist);
                    changed = true;
                end
            end
        end

        if ~changed
            break;
        end
    end
end
