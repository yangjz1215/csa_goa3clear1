function results = run_ablation_para(varargin)
    fprintf('========== 消融实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
    addParameter(p, 'n_workers', 3);
    parse(p, varargin{:});
    n_runs = p.Results.n_runs;
    map_name = p.Results.map_name;
    n_workers = p.Results.n_workers;

    [User, priorities, N_User, N_UAV, RRH, RRH_type, N_RRH, UAV_type, Ub, Lb, params, map_data] = loadExperimentContext(project_dir, map_name);

    fprintf('--- 加载固定地图: %s ---\n', map_name);
    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);
    fprintf('说明: 迭代日志跟踪当前全局标量最优解；最终汇总采用 Pareto 档案中 Utility 最大解。\n');

    variants = {
        'proposed',              'Proposed cSA-GOA (full)';
        'no_subpop',             'w/o Multi-Subpopulation';
        'no_goa_turn',           'w/o GOA Turn Operator';
        'no_goa_repulsion',      'w/o GOA Repulsion (U/V-Shape)';
        'no_pareto_leader',      'w/o Dynamic Pareto Leader';
        'no_adaptive_weight',    'w/o Adaptive Weight Rotation';
    };
    fprintf('消融实验: 变体共 %d 个\n', size(variants, 1));

    if isempty(gcp('nocreate'))
        fprintf('启动并行池 (%d workers)...\n', n_workers);
        parpool('local', n_workers);
    end

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for v_idx = 1:size(variants, 1)
        variant_name = variants{v_idx, 1};
        variant_desc = variants{v_idx, 2};
        fprintf('\n--- 运行变体: %s ---\n', variant_desc);

        scalar_utilities = zeros(n_runs, 1);
        scalar_latencies = zeros(n_runs, 1);
        scalar_energies = zeros(n_runs, 1);
        archive_utilities = zeros(n_runs, 1);
        archive_latencies = zeros(n_runs, 1);
        archive_energies = zeros(n_runs, 1);
        cov_high = zeros(n_runs, 1);
        cov_total = zeros(n_runs, 1);
        iter_counts = zeros(n_runs, 1);
        convergence_curves = cell(n_runs, 1);
        hv_values = nan(n_runs, 1);
        igd_values = nan(n_runs, 1);
        spread_values = nan(n_runs, 1);
        pareto_sizes = zeros(n_runs, 1);
        temp_pareto_fronts = cell(n_runs, 1);

        knee_utilities = zeros(n_runs, 1);
        knee_latencies = zeros(n_runs, 1);
        knee_energies = zeros(n_runs, 1);

        parfor run = 1:n_runs
            fprintf('  Worker 正在处理 Run %d/%d...\n', run, n_runs);
            stream = RandStream('mt19937ar', 'Seed', run * 1000);
            RandStream.setGlobalStream(stream);

            [~, ~, cg_curve, ~, pareto_archive, best_scalar_solution, best_utility_solution, best_knee_solution] = ...
                cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, variant_name);

            scalar_utilities(run) = best_scalar_solution.Utility;
            scalar_latencies(run) = best_scalar_solution.Latency;
            scalar_energies(run) = best_scalar_solution.Energy;

            archive_utilities(run) = best_utility_solution.Utility;
            archive_latencies(run) = best_utility_solution.Latency;
            archive_energies(run) = best_utility_solution.Energy;

            knee_utilities(run) = best_knee_solution.Utility;
            knee_latencies(run) = best_knee_solution.Latency;
            knee_energies(run) = best_knee_solution.Energy;

            bestUAV_archive = best_utility_solution.UAV_pos;
            cov_high(run) = calcCoverageWithRRH(bestUAV_archive, User(priorities >= 3, :), params.cover_radius, RRH, params.RRH_radius) * 100;
            cov_total(run) = calcCoverageWithRRH(bestUAV_archive, User, params.cover_radius, RRH, params.RRH_radius) * 100;

            iter_counts(run) = length(cg_curve);
            convergence_curves{run} = cg_curve;

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                pareto_front = archiveToFront(pareto_archive);
                pareto_sizes(run) = length(pareto_archive);
                temp_pareto_fronts{run} = pareto_front;

                sampled_front = sampleFixedSizeFront(pareto_front, 50);
                try
                    pf_norm = normalizeFront(sampled_front, priorities, N_User, N_UAV, params);
                    metrics = calculate_all_metrics(pf_norm, [], [1.1, 1.1, 1.1]);
                    hv_values(run) = metrics.hv;
                    spread_values(run) = metrics.spread;
                catch
                    hv_values(run) = NaN;
                    spread_values(run) = NaN;
                end
            else
                pareto_sizes(run) = 0;
                temp_pareto_fronts{run} = [];
            end

            fprintf('  -> Run %d 完成: ScalarUtility=%.2f, ArchiveMaxUtility=%.2f, KneeUtility=%.2f, HV=%.4f\n', ...
                run, scalar_utilities(run), archive_utilities(run), knee_utilities(run), hv_values(run));
        end

        results.(variant_name) = struct();
        results.(variant_name).description = variant_desc;
        results.(variant_name).scalar_utilities = scalar_utilities;
        results.(variant_name).scalar_latencies = scalar_latencies;
        results.(variant_name).scalar_energies = scalar_energies;
        results.(variant_name).archive_utilities = archive_utilities;
        results.(variant_name).archive_latencies = archive_latencies;
        results.(variant_name).archive_energies = archive_energies;
        results.(variant_name).knee_utilities = knee_utilities;
        results.(variant_name).knee_latencies = knee_latencies;
        results.(variant_name).knee_energies = knee_energies;
        results.(variant_name).mean_knee_utility = mean(knee_utilities);
        results.(variant_name).mean_knee_latency = mean(knee_latencies);
        results.(variant_name).mean_knee_energy = mean(knee_energies);
        results.(variant_name).cov_high = cov_high;
        results.(variant_name).cov_total = cov_total;
        results.(variant_name).iter_counts = iter_counts;
        results.(variant_name).convergence_curves = convergence_curves;
        results.(variant_name).mean_scalar_utility = mean(scalar_utilities);
        results.(variant_name).mean_scalar_latency = mean(scalar_latencies);
        results.(variant_name).mean_scalar_energy = mean(scalar_energies);
        results.(variant_name).mean_utility = mean(archive_utilities);
        results.(variant_name).std_utility = std(archive_utilities);
        results.(variant_name).mean_latency = mean(archive_latencies);
        results.(variant_name).mean_energy = mean(archive_energies);
        results.(variant_name).mean_cov_high = mean(cov_high);
        results.(variant_name).mean_cov_total = mean(cov_total);
        results.(variant_name).hv_values = hv_values;
        results.(variant_name).mean_hv = mean(hv_values, 'omitnan');
        results.(variant_name).std_hv = std(hv_values, 0, 'omitnan');
        results.(variant_name).igd_values = igd_values;
        results.(variant_name).mean_igd = mean(igd_values, 'omitnan');
        results.(variant_name).spread_values = spread_values;
        results.(variant_name).mean_spread = mean(spread_values, 'omitnan');
        results.(variant_name).pareto_fronts = temp_pareto_fronts;
        results.(variant_name).mean_pareto_size = mean(pareto_sizes);

        fprintf('  >> %s 平均结果:\n', variant_desc);
        fprintf('     日志口径(标量最优解): Utility=%.2f | Latency=%.2f s | Energy=%.2f J\n', ...
            mean(scalar_utilities), mean(scalar_latencies), mean(scalar_energies));
        fprintf('     最终汇总(档案最大Utility解): Utility=%.2f +/- %.2f\n', mean(archive_utilities), std(archive_utilities));
        fprintf('     平均时延: %.2f s\n', mean(archive_latencies));
        fprintf('     平均能耗: %.2f J\n', mean(archive_energies));
        fprintf('     平均HV: %.4f +/- %.4f\n', mean(hv_values, 'omitnan'), std(hv_values, 0, 'omitnan'));
    end

    results = finalizeNormalizedMetrics(results, priorities, N_User, N_UAV, params, n_runs, '消融实验');

    results_file = fullfile(project_dir, 'experiments', ['ablation_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results', 'map_data', 'params');
    fprintf('\n结果已保存至: %s\n', results_file);
end

function [User, priorities, N_User, N_UAV, RRH, RRH_type, N_RRH, UAV_type, Ub, Lb, params, map_data] = loadExperimentContext(project_dir, map_name)
    maps_dir = fullfile(project_dir, 'maps');
    map_file = fullfile(maps_dir, [map_name, '.mat']);
    if ~exist(map_file, 'file')
        error('地图文件不存在: %s', map_file);
    end

    map_data = load(map_file);
    User = map_data.User;
    priorities = map_data.priorities;
    N_User = map_data.N_User;

    if size(User, 1) == 200
        N_UAV = 8;
        config_suffix = 'Small';
    elseif size(User, 1) == 500
        N_UAV = 15;
        config_suffix = 'Medium';
    else
        N_UAV = 25;
        config_suffix = 'Large';
    end

    config_file = fullfile(maps_dir, ['Map_', config_suffix, '_Config.mat']);
    if exist(config_file, 'file')
        config_data = load(config_file);
        RRH = config_data.RRH;
        RRH_type = config_data.RRH_type;
        N_RRH = config_data.N_RRH;
        UAV_type = config_data.UAV_type;
        params.D = config_data.D;
        params.C = config_data.C;
        params.DT = config_data.DT;
    else
        N_RRH = 10;
        RRH = rand(N_RRH, 2) * 1000;
        RRH_type = zeros(N_RRH, 1);
        UAV_type = zeros(N_UAV, 1);
        params.D = [];
        params.C = [];
        params.DT = [];
    end

    Lb = [0, 0];
    Ub = [1000, 1000];
    params.RRH = RRH;
    params.RRH_type = RRH_type;
    params.UAV_type = UAV_type;
    params.max_latency = 1.0;
    params.B_total = 10e6;
    params.enable_bilevel = true;
    params.F_total = 10e9;
    params.noise = 1e-13;
    params.P_tx = 0.5;
    params.PtxU = 0.5;
    params.B_total_relay = 20e6;
    params.f_BBU = 50e9;
    params.kappa = 1e-28;
    params.k_move = 10;
    params.cover_radius = 150;
    params.RRH_radius = 150;
    params.E_max = 50000;
    params.energy_norm_max = 15000;
    params.D_UU = 10;
    params.D_RU = 10;
    params.FES_max = 300;
    params.K = 40;
    params.G_weights = [0.4, 0.3, 0.3];
    params.subpop_params = struct();
    params.subpop_params.mu0 = [500, 500];
    params.subpop_params.sigma0 = [150, 150; 120, 120; 80, 80];
    params.subpop_params.sigma_min = [5, 8, 3];
    params.subpop_params.q = [0.6, 0.5, 0.4];
    params.subpop_params.beta = [0.8, 0.7, 0.6];

    if ~isfield(params, 'D') || isempty(params.D) || length(params.D) ~= N_User
        params.D = ones(N_User, 1) * 2e6;
    end
    if ~isfield(params, 'C') || isempty(params.C) || length(params.C) ~= N_User
        params.C = ones(N_User, 1) * 0.5e9;
    end

    hyperparam_file = fullfile(project_dir, 'experiments', 'best_algo_hyperparams.mat');
    if exist(hyperparam_file, 'file')
        opt_params = load(hyperparam_file);
        params.K = opt_params.best_K;
        params.subpop_params.q = [opt_params.best_q, opt_params.best_q, opt_params.best_q];
        fprintf('[动态注入] 已加载最优超参数 K=%d, q=%.1f\n', params.K, opt_params.best_q);
        if isfield(opt_params, 'best_beta')
            params.subpop_params.beta = [opt_params.best_beta, opt_params.best_beta, opt_params.best_beta];
            fprintf('[动态注入] beta=%.2f\n', opt_params.best_beta);
        end
    end
end

function pareto_front = archiveToFront(pareto_archive)
    pareto_front = zeros(length(pareto_archive), 3);
    for p_idx = 1:length(pareto_archive)
        pareto_front(p_idx, 1) = pareto_archive(p_idx).Utility;
        pareto_front(p_idx, 2) = pareto_archive(p_idx).Latency;
        pareto_front(p_idx, 3) = pareto_archive(p_idx).Energy;
    end
end

function pf_norm = normalizeFront(pareto_front, priorities, N_User, N_UAV, params)
    max_U = sum(priorities);
    max_L = N_User * params.max_latency;
    max_E = params.E_max * N_UAV;
    pf_norm = zeros(size(pareto_front));
    pf_norm(:, 1) = (max_U - pareto_front(:, 1)) / (max_U + 1e-6);
    pf_norm(:, 2) = pareto_front(:, 2) / (max_L + 1e-6);
    pf_norm(:, 3) = pareto_front(:, 3) / (max_E + 1e-6);
end

function results = finalizeNormalizedMetrics(results, priorities, N_User, N_UAV, params, n_runs, title_text)
    fprintf('\n========== 开始执行 %s 归一化 + Leave-One-Out IGD ==========\n', title_text);

    field_names = fieldnames(results);
    valid_names = {};

    % Step 1: 归一化所有前沿，按算法收集
    all_fronts = struct();
    for i = 1:length(field_names)
        name = field_names{i};
        if isstruct(results.(name)) && isfield(results.(name), 'pareto_fronts')
            valid_names{end + 1} = name; %#ok<AGROW>
            all_fronts.(name) = [];
            for run = 1:length(results.(name).pareto_fronts)
                pf = results.(name).pareto_fronts{run};
                if ~isempty(pf)
                    pf_norm = normalizeFront(pf, priorities, N_User, N_UAV, params);
                    all_fronts.(name) = [all_fronts.(name); pf_norm]; %#ok<AGROW>
                    results.(name).pareto_fronts_norm{run} = pf_norm;
                else
                    results.(name).pareto_fronts_norm{run} = [];
                end
            end
        end
    end

    ref_point_norm = [1.1, 1.1, 1.1];

    fprintf('\n========== %s最终汇总表格 (统一口径) ==========\n', title_text);
    fprintf('%-14s | %-10s | %-12s | %-12s | %-13s | %-13s | %-8s\n', ...
        'Name', 'Utility', 'Latency(s)', 'Energy(J)', 'HV(Norm)↑', 'IGD(Norm)↓', 'Spread↑');
    fprintf('%s\n', repmat('-', 1, 120));

    % Step 2: Leave-one-out 计算每个算法的指标
    for i = 1:length(valid_names)
        name = valid_names{i};
        r = results.(name);

        % 构建排除自身的参考前沿
        other_points = [];
        for j = 1:length(valid_names)
            if j ~= i
                other_points = [other_points; all_fronts.(valid_names{j})]; %#ok<AGROW>
            end
        end
        if ~isempty(other_points)
            ref_front = extractNonDominated(unique(other_points, 'rows'));
        else
            ref_front = [];
        end

        hvs = nan(n_runs, 1);
        igds = nan(n_runs, 1);
        gds = nan(n_runs, 1);
        spreads = nan(n_runs, 1);

        for run = 1:length(r.pareto_fronts_norm)
            pf_norm = r.pareto_fronts_norm{run};
            if ~isempty(pf_norm)
                metrics = calculate_all_metrics(pf_norm, ref_front, ref_point_norm);
                hvs(run) = metrics.hv;
                igds(run) = metrics.igd;
                gds(run) = metrics.gd;
                spreads(run) = metrics.spread;
            end
        end

        results.(name).mean_hv_norm = mean(hvs, 'omitnan');
        results.(name).std_hv_norm = std(hvs, 0, 'omitnan');
        results.(name).mean_igd_norm = mean(igds, 'omitnan');
        results.(name).std_igd_norm = std(igds, 0, 'omitnan');
        results.(name).mean_gd_norm = mean(gds, 'omitnan');
        results.(name).std_gd_norm = std(gds, 0, 'omitnan');
        results.(name).mean_spread_norm = mean(spreads, 'omitnan');
        results.(name).std_spread_norm = std(spreads, 0, 'omitnan');

        results.(name).hv_values = hvs;
        results.(name).igd_values = igds;
        results.(name).gd_values = gds;
        results.(name).spread_values = spreads;

        fprintf('%-14s | %-10.2f | %-12.2f | %-12.2f | HV=%.4f±%.4f | IGD=%.4f±%.4f | GD=%.4f±%.4f | Spread=%.4f\n', ...
            name, r.mean_utility, r.mean_latency, r.mean_energy, mean(hvs,'omitnan'), std(hvs,0,'omitnan'), ...
            mean(igds,'omitnan'), std(igds,0,'omitnan'), mean(gds,'omitnan'), std(gds,0,'omitnan'), mean(spreads,'omitnan'));
    end
    fprintf('%s\n', repmat('-', 1, 120));
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
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

function pf = extractNonDominated(points)
    if isempty(points)
        pf = [];
        return;
    end

    n = size(points, 1);
    is_dominated = false(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j && all(points(j, :) <= points(i, :)) && any(points(j, :) < points(i, :))
                is_dominated(i) = true;
                break;
            end
        end
    end
    pf = points(~is_dominated, :);
end

function sampled_front = sampleFixedSizeFront(pareto_front, target_size)
    if isempty(pareto_front)
        sampled_front = [];
        return;
    end

    nd_front = extractNonDominated(pareto_front);

    if size(nd_front, 1) <= target_size
        sampled_front = nd_front;
        return;
    end

    n = size(nd_front, 1);
    n_obj = size(nd_front, 2);
    crowding = zeros(n, 1);

    for m = 1:n_obj
        [sorted_obj, sort_idx] = sort(nd_front(:, m));
        crowding(sort_idx(1)) = inf;
        crowding(sort_idx(end)) = inf;

        f_min = sorted_obj(1);
        f_max = sorted_obj(end);

        if f_max == f_min
            continue;
        end

        for i = 2:(n-1)
            crowding(sort_idx(i)) = crowding(sort_idx(i)) + ...
                (sorted_obj(i+1) - sorted_obj(i-1)) / (f_max - f_min);
        end
    end

    [~, sorted_crowd_idx] = sort(crowding, 'descend');
    selected_idx = sorted_crowd_idx(1:target_size);
    sampled_front = nd_front(selected_idx, :);
end
