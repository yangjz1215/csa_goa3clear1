function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = NSGA2_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
    end

    if isfield(params, 'K')
        pop_size = params.K;
    else
        pop_size = 40;
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
    n_obj = 3;

    pc = 0.9;
    pm = 1.0 / n_vars;
    eta_c = 20;
    eta_m = 20;

    population = zeros(pop_size, n_vars);
    objectives = zeros(pop_size, n_obj);

    center_point = [500, 500];
    jitter = 10;

    for i = 1:pop_size
        for j = 1:N_UAV
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            population(i, (j-1)*2+1) = max(Lb(1), min(Ub(1), init_x));
            population(i, (j-1)*2+2) = max(Lb(2), min(Ub(2), init_y));
        end
    end

    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

    for i = 1:pop_size
        uav_pos = reshape(population(i, :), N_UAV, 2);
        if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
            uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
            population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
        end
        population(i, :) = max(0, min(1000, population(i, :)));

        [util, lat, nrg, ~] = calcMEC_Objectives(uav_pos, User, priorities, params);
        objectives(i, 1) = -util;
        objectives(i, 2) = lat;
        objectives(i, 3) = nrg;

        pareto_archive = updateParetoArchive3D(pareto_archive, uav_pos, util, lat, nrg);
    end

    cg_curve = zeros(1, max_iter);
    best_util_so_far = -inf;
    bestUAV = zeros(N_UAV, 2);
    best_energy = inf;

    for i = 1:pop_size
        if -objectives(i, 1) > best_util_so_far
            best_util_so_far = -objectives(i, 1);
            bestUAV = reshape(population(i, :), N_UAV, 2);
            best_energy = objectives(i, 3);
        end
    end
    cg_curve(1) = best_util_so_far;

    [fronts, rank] = fastNonDominatedSort(objectives);
    crowd = crowdingDistance(objectives, fronts);

    for iter = 2:max_iter
        offspring = zeros(pop_size, n_vars);

        for i = 1:2:pop_size
            p1_idx = binaryTournament(rank, crowd);
            p2_idx = binaryTournament(rank, crowd);

            parent1 = population(p1_idx, :);
            parent2 = population(p2_idx, :);

            if rand() < pc
                [child1, child2] = sbxCrossover(parent1, parent2, Lb, Ub, eta_c);
            else
                child1 = parent1;
                child2 = parent2;
            end

            child1 = polynomialMutation(child1, Lb, Ub, pm, eta_m);
            child2 = polynomialMutation(child2, Lb, Ub, pm, eta_m);

            offspring(i, :) = child1;
            if i + 1 <= pop_size
                offspring(i+1, :) = child2;
            end
        end

        for i = 1:pop_size
            uav_pos = reshape(offspring(i, :), N_UAV, 2);
            if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
                uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
                offspring(i, :) = reshape(uav_pos, 1, N_UAV * 2);
            end
            offspring(i, :) = max(0, min(1000, offspring(i, :)));
        end

        off_obj = zeros(pop_size, n_obj);
        for i = 1:pop_size
            uav_pos = reshape(offspring(i, :), N_UAV, 2);
            [util, lat, nrg, ~] = calcMEC_Objectives(uav_pos, User, priorities, params);
            off_obj(i, 1) = -util;
            off_obj(i, 2) = lat;
            off_obj(i, 3) = nrg;

            pareto_archive = updateParetoArchive3D(pareto_archive, uav_pos, util, lat, nrg);

            if util > best_util_so_far
                best_util_so_far = util;
                bestUAV = uav_pos;
                best_energy = nrg;
            end
        end

        combined_pop = [population; offspring];
        combined_obj = [objectives; off_obj];

        [population, objectives] = nsga2Selection(combined_pop, combined_obj, pop_size);
        [fronts, rank] = fastNonDominatedSort(objectives);
        crowd = crowdingDistance(objectives, fronts);

        cg_curve(iter) = best_util_so_far;

        if mod(iter, 50) == 0
            fprintf('NSGA-II iter %d/%d, Best utility: %.2f, Archive size: %d\n', ...
                iter, max_iter, best_util_so_far, length(pareto_archive));
        end
    end

    if ~isempty(pareto_archive) && length(pareto_archive) >= 1
        arch_utils = [pareto_archive.Utility];
        [best_fit, max_u_idx] = max(arch_utils);
        bestUAV = pareto_archive(max_u_idx).UAV_pos;
        best_energy = pareto_archive(max_u_idx).Energy;
    else
        best_fit = best_util_so_far;
    end
end

function idx = binaryTournament(rank, crowd)
    n = length(rank);
    i1 = randi(n);
    i2 = randi(n);

    if rank(i1) < rank(i2)
        idx = i1;
    elseif rank(i2) < rank(i1)
        idx = i2;
    elseif crowd(i1) > crowd(i2)
        idx = i1;
    else
        idx = i2;
    end
end

function [fronts, rank] = fastNonDominatedSort(objectives)
    n = size(objectives, 1);

    dominated_count = zeros(n, 1);
    dominates_list = cell(n, 1);

    for i = 1:n
        for j = i+1:n
            if all(objectives(i, :) <= objectives(j, :)) && any(objectives(i, :) < objectives(j, :))
                dominates_list{i} = [dominates_list{i}, j];
                dominated_count(j) = dominated_count(j) + 1;
            elseif all(objectives(j, :) <= objectives(i, :)) && any(objectives(j, :) < objectives(i, :))
                dominates_list{j} = [dominates_list{j}, i];
                dominated_count(i) = dominated_count(i) + 1;
            end
        end
    end

    fronts = {};
    current_front = find(dominated_count == 0)';
    rank = zeros(n, 1);

    while ~isempty(current_front)
        fronts{end+1} = current_front;
        rank(current_front) = length(fronts);
        next_front = [];
        for i = current_front
            for j = dominates_list{i}
                dominated_count(j) = dominated_count(j) - 1;
                if dominated_count(j) == 0
                    next_front = [next_front, j];
                end
            end
        end
        current_front = next_front;
    end
end

function crowd = crowdingDistance(objectives, fronts)
    n = size(objectives, 1);
    m = size(objectives, 2);
    crowd = zeros(n, 1);

    for f = 1:length(fronts)
        front = fronts{f};
        nf = length(front);

        if nf < 3
            crowd(front) = inf;
            continue;
        end

        for obj = 1:m
            [~, sorted_idx] = sort(objectives(front, obj));
            sorted_front = front(sorted_idx);

            obj_range = objectives(sorted_front(end), obj) - objectives(sorted_front(1), obj);
            if obj_range == 0, obj_range = 1e-6; end

            crowd(sorted_front(1)) = inf;
            crowd(sorted_front(end)) = inf;

            for k = 2:nf-1
                crowd(sorted_front(k)) = crowd(sorted_front(k)) + ...
                    (objectives(sorted_front(k+1), obj) - objectives(sorted_front(k-1), obj)) / obj_range;
            end
        end
    end
end

function [new_pop, new_obj] = nsga2Selection(combined_pop, combined_obj, pop_size)
    [fronts, rank] = fastNonDominatedSort(combined_obj);
    crowd = crowdingDistance(combined_obj, fronts);

    new_pop = zeros(pop_size, size(combined_pop, 2));
    new_obj = zeros(pop_size, size(combined_obj, 2));
    count = 0;
    f = 1;

    while f <= length(fronts) && count + length(fronts{f}) <= pop_size
        front = fronts{f};
        new_pop(count+1:count+length(front), :) = combined_pop(front, :);
        new_obj(count+1:count+length(front), :) = combined_obj(front, :);
        count = count + length(front);
        f = f + 1;
    end

    if count < pop_size && f <= length(fronts)
        last_front = fronts{f};
        [~, sorted_idx] = sort(crowd(last_front), 'descend');
        remaining = pop_size - count;
        selected = last_front(sorted_idx(1:remaining));
        new_pop(count+1:end, :) = combined_pop(selected, :);
        new_obj(count+1:end, :) = combined_obj(selected, :);
    end
end

function [child1, child2] = sbxCrossover(parent1, parent2, Lb, Ub, eta_c)
    n = length(parent1);
    child1 = zeros(1, n);
    child2 = zeros(1, n);

    for i = 1:n
        if rand() <= 0.5
            if abs(parent2(i) - parent1(i)) > 1e-14
                if parent1(i) < parent2(i)
                    y1 = parent1(i);
                    y2 = parent2(i);
                else
                    y1 = parent2(i);
                    y2 = parent1(i);
                end

                yl = Lb(mod(i-1, 2) + 1);
                yu = Ub(mod(i-1, 2) + 1);

                beta = 1.0 + (2.0 * (y1 - yl) / max(y2 - y1, 1e-14));
                alpha = 2.0 - beta^(-(eta_c + 1.0));
                u = rand();
                if u <= 1.0 / alpha
                    betaq = (u * alpha)^(1.0 / (eta_c + 1.0));
                else
                    betaq = (1.0 / (2.0 - u * alpha))^(1.0 / (eta_c + 1.0));
                end
                c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

                beta = 1.0 + (2.0 * (yu - y2) / max(y2 - y1, 1e-14));
                alpha = 2.0 - beta^(-(eta_c + 1.0));
                u = rand();
                if u <= 1.0 / alpha
                    betaq = (u * alpha)^(1.0 / (eta_c + 1.0));
                else
                    betaq = (1.0 / (2.0 - u * alpha))^(1.0 / (eta_c + 1.0));
                end
                c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                c1 = max(yl, min(yu, c1));
                c2 = max(yl, min(yu, c2));

                if rand() <= 0.5
                    child1(i) = c2;
                    child2(i) = c1;
                else
                    child1(i) = c1;
                    child2(i) = c2;
                end
            else
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        else
            child1(i) = parent1(i);
            child2(i) = parent2(i);
        end
    end
end

function child = polynomialMutation(parent, Lb, Ub, pm, eta_m)
    n = length(parent);
    child = parent;

    for i = 1:n
        if rand() < pm
            y = parent(i);
            yl = Lb(mod(i-1, 2) + 1);
            yu = Ub(mod(i-1, 2) + 1);

            delta1 = (y - yl) / max(yu - yl, 1e-14);
            delta2 = (yu - y) / max(yu - yl, 1e-14);

            u = rand();
            if u <= 0.5
                xy = 1.0 - delta1;
                val = 2.0 * u + (1.0 - 2.0 * u) * xy^(eta_m + 1.0);
                deltaq = val^(1.0 / (eta_m + 1.0)) - 1.0;
            else
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - u) + 2.0 * (u - 0.5) * xy^(eta_m + 1.0);
                deltaq = 1.0 - val^(1.0 / (eta_m + 1.0));
            end

            y = y + deltaq * (yu - yl);
            child(i) = max(yl, min(yu, y));
        end
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