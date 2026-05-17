function results = run_comparison_para(varargin)
    fprintf('========== 对比实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'comparison_algorithms')));
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
    params.enable_early_stop = false;
    params.enable_smart_stop = false;
    params.enable_migration_log = false;

    fprintf('--- 加载固定地图: %s ---\n', map_name);
    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);
    fprintf('说明: cSA-GOA 迭代日志跟踪当前全局标量最优解；最终汇总统一采用各算法档案中的最大 Utility 解。\n');

    algorithms = {
        'cSA_GOA', 'cSA-GOA (Proposed)';
        'PSO', 'PSO (Particle Swarm Optimization)';
        'GA', 'GA (Genetic Algorithm)';
        'GOA', 'GOA (Grasshopper Optimization)';
        'cSA', 'cSA (Compact Sine Algorithm)';
        'NSGA2', 'NSGA-II (Non-dominated Sorting GA)'
    };

    if isempty(gcp('nocreate'))
        fprintf('启动并行池 (%d workers)...\n', n_workers);
        parpool('local', n_workers);
    end

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for alg_idx = 1:size(algorithms, 1)
        alg_name = algorithms{alg_idx, 1};
        alg_desc = algorithms{alg_idx, 2};
        fprintf('\n--- 运行算法: %s (并行) ---\n', alg_desc);

        scalar_utilities = nan(n_runs, 1);
        scalar_latencies = nan(n_runs, 1);
        scalar_energies = nan(n_runs, 1);
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
        pareto_fronts_cell = cell(n_runs, 1);
        runtimes = zeros(n_runs, 1);
        success_rates = zeros(n_runs, 1);

        parfor run = 1:n_runs
            stream = RandStream('mt19937ar', 'Seed', run * 200 + alg_idx * 1000);
            RandStream.setGlobalStream(stream);
            t_start = tic;

            bestUAV = zeros(N_UAV, 2);
            cg_curve = zeros(1, params.FES_max);
            pareto_archive = [];
            best_scalar_solution = [];
            best_utility_solution = [];

            switch alg_name
                case 'cSA_GOA'
                    [~, bestUAV, cg_curve, ~, ~, ~, ~, ~, ~, pareto_archive, best_scalar_solution, best_utility_solution] = ...
                        cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'PSO'
                    [~, bestUAV, cg_curve, ~, pareto_archive] = ...
                        PSO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GA'
                    [~, bestUAV, cg_curve, ~, pareto_archive] = ...
                        GA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GOA'
                    [~, bestUAV, cg_curve, ~, pareto_archive] = ...
                        GOA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'cSA'
                    [~, bestUAV, cg_curve, ~, pareto_archive] = ...
                        cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'NSGA2'
                    [~, bestUAV, cg_curve, ~, pareto_archive] = ...
                        NSGA2_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
            end

            runtimes(run) = toc(t_start);

            if isempty(best_scalar_solution)
                best_scalar_solution = buildSolutionFromUAV(bestUAV, User, priorities, params, 'scalar_fallback');
            end
            if isempty(best_utility_solution)
                best_utility_solution = selectUtilitySolutionFromArchive(pareto_archive, bestUAV, User, priorities, params);
            end

            scalar_utilities(run) = best_scalar_solution.Utility;
            scalar_latencies(run) = best_scalar_solution.Latency;
            scalar_energies(run) = best_scalar_solution.Energy;
            archive_utilities(run) = best_utility_solution.Utility;
            archive_latencies(run) = best_utility_solution.Latency;
            archive_energies(run) = best_utility_solution.Energy;
            success_rates(run) = best_utility_solution.SuccessRate;

            bestUAV_archive = best_utility_solution.UAV_pos;
            cov_high(run) = calcCoverageWithRRH(bestUAV_archive, User(priorities >= 3, :), params.cover_radius, RRH, params.RRH_radius) * 100;
            cov_total(run) = calcCoverageWithRRH(bestUAV_archive, User, params.cover_radius, RRH, params.RRH_radius) * 100;

            cg_curve = sanitizeConvergenceCurve(cg_curve);
            iter_counts(run) = length(cg_curve);
            convergence_curves{run} = cg_curve;

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                pareto_front = archiveToFront(pareto_archive);
                pareto_sizes(run) = length(pareto_archive);
                pareto_fronts_cell{run} = pareto_front;

                try
                    pf_norm = normalizeFront(pareto_front, priorities, N_User, N_UAV, params);
                    metrics = calculate_all_metrics(pf_norm, [], [1.1, 1.1, 1.1]);
                    hv_values(run) = metrics.hv;
                    spread_values(run) = metrics.spread;
                catch
                    hv_values(run) = NaN;
                    spread_values(run) = NaN;
                end
            else
                pareto_fronts_cell{run} = [];
                pareto_sizes(run) = 0;
            end

            fprintf('  -> Run %d 完成: ScalarUtility=%.2f, ArchiveMaxUtility=%.2f, HV=%.4f\n', ...
                run, scalar_utilities(run), archive_utilities(run), hv_values(run));
        end

        results.(alg_name) = struct();
        results.(alg_name).description = alg_desc;
        results.(alg_name).scalar_utilities = scalar_utilities;
        results.(alg_name).scalar_latencies = scalar_latencies;
        results.(alg_name).scalar_energies = scalar_energies;
        results.(alg_name).archive_utilities = archive_utilities;
        results.(alg_name).archive_latencies = archive_latencies;
        results.(alg_name).archive_energies = archive_energies;
        results.(alg_name).cov_high = cov_high;
        results.(alg_name).cov_total = cov_total;
        results.(alg_name).iter_counts = iter_counts;
        results.(alg_name).convergence_curves = convergence_curves;
        results.(alg_name).convergence_metric = 'best_so_far_utility';
        results.(alg_name).convergence_unit = 'priority_sum';
        results.(alg_name).max_generations = params.FES_max;
        results.(alg_name).mean_scalar_utility = mean(scalar_utilities, 'omitnan');
        results.(alg_name).mean_scalar_latency = mean(scalar_latencies, 'omitnan');
        results.(alg_name).mean_scalar_energy = mean(scalar_energies, 'omitnan');
        results.(alg_name).mean_utility = mean(archive_utilities);
        results.(alg_name).std_utility = std(archive_utilities);
        results.(alg_name).mean_latency = mean(archive_latencies);
        results.(alg_name).mean_energy = mean(archive_energies);
        results.(alg_name).mean_cov_high = mean(cov_high);
        results.(alg_name).mean_cov_total = mean(cov_total);
        results.(alg_name).hv_values = hv_values;
        results.(alg_name).mean_hv = mean(hv_values, 'omitnan');
        results.(alg_name).std_hv = std(hv_values, 0, 'omitnan');
        results.(alg_name).igd_values = igd_values;
        results.(alg_name).mean_igd = mean(igd_values, 'omitnan');
        results.(alg_name).spread_values = spread_values;
        results.(alg_name).mean_spread = mean(spread_values, 'omitnan');
        results.(alg_name).mean_pareto_size = mean(pareto_sizes);
        results.(alg_name).pareto_fronts = pareto_fronts_cell;
        results.(alg_name).runtimes = runtimes;
        results.(alg_name).mean_runtime = mean(runtimes);
        results.(alg_name).success_rates = success_rates;
        results.(alg_name).mean_success_rate = mean(success_rates);

        fprintf('  >> %s 平均结果:\n', alg_desc);
        fprintf('     日志口径(标量最优解): Utility=%.2f | Latency=%.2f s | Energy=%.2f J\n', ...
            mean(scalar_utilities, 'omitnan'), mean(scalar_latencies, 'omitnan'), mean(scalar_energies, 'omitnan'));
        fprintf('     最终汇总(档案最大Utility解): Utility=%.2f +/- %.2f\n', mean(archive_utilities), std(archive_utilities));
        fprintf('     平均时延: %.2f s\n', mean(archive_latencies));
        fprintf('     平均能耗: %.2f J\n', mean(archive_energies));
        fprintf('     平均HV: %.4f +/- %.4f\n', mean(hv_values, 'omitnan'), std(hv_values, 0, 'omitnan'));
        fprintf('     平均Runtime: %.2f s\n', mean(runtimes));
    end

    results = finalizeNormalizedMetrics(results, priorities, N_User, N_UAV, params, n_runs, '对比实验');

    fprintf('\n========== Wilcoxon Rank Sum Test (cSA-GOA vs. Others) ==========\n');
    wilcoxon_table_hv = computeWilcoxonTable(results, algorithms, 'hv_values');
    wilcoxon_table_igd = computeWilcoxonTable(results, algorithms, 'igd_values');
    wilcoxon_table_spread = computeWilcoxonTable(results, algorithms, 'spread_values');
    results.wilcoxon_table_hv = wilcoxon_table_hv;
    results.wilcoxon_table_igd = wilcoxon_table_igd;
    results.wilcoxon_table_spread = wilcoxon_table_spread;

    results_file = fullfile(project_dir, 'experiments', ['comparison_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results', 'map_data', 'params');
    fprintf('\n结果已保存至: %s\n', results_file);
end

function curve = sanitizeConvergenceCurve(curve)
    curve = double(curve(:)');
    if isempty(curve)
        return;
    end

    curve(~isfinite(curve)) = NaN;
    first_valid = find(~isnan(curve), 1, 'first');
    if isempty(first_valid)
        curve = zeros(size(curve));
        return;
    end

    curve(1:first_valid-1) = curve(first_valid);
    for idx = first_valid + 1:length(curve)
        if isnan(curve(idx))
            curve(idx) = curve(idx - 1);
        end
    end

    curve = cummax(curve);
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
        fprintf('🎯 [动态注入] 已加载最优超参数 K=%d, q=%.1f\n', params.K, opt_params.best_q);
        if isfield(opt_params, 'best_phase_w_progress')
            params.phase_w_progress = opt_params.best_phase_w_progress;
            params.phase_w_cov = opt_params.best_phase_w_cov;
            params.phase_w_inner = opt_params.best_phase_w_inner;
            fprintf('🎯 [动态注入] phi_t 权重: progress=%.3f, cov=%.3f, inner=%.3f\n', ...
                opt_params.best_phase_w_progress, opt_params.best_phase_w_cov, opt_params.best_phase_w_inner);
        end
        if isfield(opt_params, 'best_beta')
            params.subpop_params.beta = [opt_params.best_beta, opt_params.best_beta, opt_params.best_beta];
            fprintf('🎯 [动态注入] beta=%.2f\n', opt_params.best_beta);
        end
    end
end

function solution = buildSolutionFromUAV(bestUAV, User, priorities, params, label)
    [util, lat, nrg, success_rate] = calcMEC_Objectives(bestUAV, User, priorities, params);
    solution = struct('label', label, 'UAV_pos', bestUAV, 'Utility', util, 'Latency', lat, 'Energy', nrg, 'SuccessRate', success_rate);
end

function solution = selectUtilitySolutionFromArchive(pareto_archive, fallback_uav, User, priorities, params)
    if ~isempty(pareto_archive)
        arch_util = [pareto_archive.Utility];
        [~, max_idx] = max(arch_util);
        solution = buildSolutionFromUAV(pareto_archive(max_idx).UAV_pos, User, priorities, params, 'archive_max_utility');
        solution.Utility = pareto_archive(max_idx).Utility;
        solution.Latency = pareto_archive(max_idx).Latency;
        solution.Energy = pareto_archive(max_idx).Energy;
    else
        solution = buildSolutionFromUAV(fallback_uav, User, priorities, params, 'archive_max_utility_fallback');
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
        spreads = nan(n_runs, 1);

        for run = 1:length(r.pareto_fronts_norm)
            pf_norm = r.pareto_fronts_norm{run};
            if ~isempty(pf_norm)
                metrics = calculate_all_metrics(pf_norm, ref_front, ref_point_norm);
                hvs(run) = metrics.hv;
                igds(run) = metrics.igd;
                spreads(run) = metrics.spread;
            end
        end

        results.(name).mean_hv_norm = mean(hvs, 'omitnan');
        results.(name).std_hv_norm = std(hvs, 0, 'omitnan');
        results.(name).mean_igd_norm = mean(igds, 'omitnan');
        results.(name).std_igd_norm = std(igds, 0, 'omitnan');
        results.(name).mean_spread_norm = mean(spreads, 'omitnan');
        results.(name).std_spread_norm = std(spreads, 0, 'omitnan');

        results.(name).igd_values = igds;
        results.(name).spread_values = spreads;
        results.(name).hv_values = hvs;

        fprintf('%-14s | %-10.2f | %-12.2f | %-12.2f | %-5.4f±%-5.4f | %-5.4f±%-5.4f | %-8.4f\n', ...
            name, r.mean_utility, r.mean_latency, r.mean_energy, mean(hvs, 'omitnan'), std(hvs, 0, 'omitnan'), mean(igds, 'omitnan'), std(igds, 0, 'omitnan'), mean(spreads, 'omitnan'));
    end
    fprintf('%s\n', repmat('-', 1, 120));
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos, 1)
        user = User_pos(i, :);
        dist_UAV = min(sqrt(sum((UAV_pos - repmat(user, size(UAV_pos, 1), 1)).^2, 2)));
        dist_RRH = min(sqrt(sum((RRH - repmat(user, size(RRH, 1), 1)).^2, 2)));
        if dist_UAV <= UAV_radius || dist_RRH <= RRH_radius
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

function wilcoxon_table = computeWilcoxonTable(results, algorithms, metric_field)
    if nargin < 3
        metric_field = 'hv_values';
    end

    n_algs = size(algorithms, 1);
    alg_names = algorithms(:, 1);
    alg_labels = algorithms(:, 2);

    ref_name = 'cSA_GOA';
    if ~isfield(results, ref_name) || ~isfield(results.(ref_name), metric_field)
        warning('cSA_GOA %s data not found, skipping Wilcoxon test.', metric_field);
        wilcoxon_table = [];
        return;
    end

    ref_vals = results.(ref_name).(metric_field)(:);
    ref_vals = ref_vals(~isnan(ref_vals));

    fprintf('\n--- Wilcoxon on %s ---\n', metric_field);
    fprintf('%-40s | %-12s | %-12s | %-14s\n', 'Comparison', 'p-value', 'h (p<0.05)', 'Significant');
    fprintf('%s\n', repmat('-', 1, 85));

    wilcoxon_table = cell(n_algs, 5);
    row = 0;

    for a = 1:n_algs
        alg_name = alg_names{a};
        if strcmp(alg_name, ref_name), continue; end
        if ~isfield(results, alg_name) || ~isfield(results.(alg_name), metric_field), continue; end

        cmp_vals = results.(alg_name).(metric_field)(:);
        cmp_vals = cmp_vals(~isnan(cmp_vals));

        if length(ref_vals) < 3 || length(cmp_vals) < 3
            p_val = NaN;
            h = NaN;
        else
            [p_val, h] = wilcoxon_test(ref_vals, cmp_vals);
        end

        sig_str = '';
        if ~isnan(h) && h == 1
            sig_str = 'Yes (p<0.05)';
        elseif ~isnan(h)
            sig_str = 'No (p>=0.05)';
        else
            sig_str = 'N/A';
        end

        row = row + 1;
        wilcoxon_table{row, 1} = ['cSA-GOA vs. ', alg_labels{a}];
        wilcoxon_table{row, 2} = p_val;
        wilcoxon_table{row, 3} = h;
        wilcoxon_table{row, 4} = sig_str;
        wilcoxon_table{row, 5} = alg_name;

        fprintf('%-40s | %-12.2e | %-12d | %-14s\n', ...
            wilcoxon_table{row, 1}, p_val, h, sig_str);
    end

    wilcoxon_table = wilcoxon_table(1:row, :);
    fprintf('%s\n', repmat('-', 1, 85));
end
