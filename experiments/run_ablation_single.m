function result = run_ablation_single(variant_name, n_runs, map_name)
% run_ablation_single - 运行单个消融变体，结果追加保存到临时文件
% 用法: run_ablation_single('no_subpop', 5, 'Map1_Medium')
%       run_ablation_single('no_goa_turn', 5, 'Map1_Medium')

    fprintf('========== 单变体消融: %s (%d次, %s) ==========\n', variant_name, n_runs, map_name);

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    [User, priorities, N_User, N_UAV, RRH, RRH_type, N_RRH, UAV_type, Ub, Lb, params, map_data] = loadExperimentContext(project_dir, map_name);

    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);

    scalar_utilities = zeros(n_runs, 1);
    scalar_latencies = zeros(n_runs, 1);
    scalar_energies = zeros(n_runs, 1);
    archive_utilities = zeros(n_runs, 1);
    archive_latencies = zeros(n_runs, 1);
    archive_energies = zeros(n_runs, 1);
    knee_utilities = zeros(n_runs, 1);
    knee_latencies = zeros(n_runs, 1);
    knee_energies = zeros(n_runs, 1);
    cov_high = zeros(n_runs, 1);
    cov_total = zeros(n_runs, 1);
    convergence_curves = cell(n_runs, 1);
    hv_values = nan(n_runs, 1);
    pareto_sizes = zeros(n_runs, 1);
    temp_pareto_fronts = cell(n_runs, 1);

    for run = 1:n_runs
        fprintf('  Run %d/%d...\n', run, n_runs);
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

        convergence_curves{run} = cg_curve;

        if ~isempty(pareto_archive) && length(pareto_archive) > 1
            pareto_front = archiveToFront(pareto_archive);
            pareto_sizes(run) = length(pareto_archive);
            temp_pareto_fronts{run} = pareto_front;
        else
            pareto_sizes(run) = 0;
            temp_pareto_fronts{run} = [];
        end

        fprintf('  -> Run %d: ArchiveUtility=%.2f, HV=%.4f\n', run, archive_utilities(run), hv_values(run));
    end

    result = struct();
    result.variant = variant_name;
    result.mean_utility = mean(archive_utilities);
    result.std_utility = std(archive_utilities);
    result.mean_latency = mean(archive_latencies);
    result.mean_energy = mean(archive_energies);
    result.mean_knee_utility = mean(knee_utilities);
    result.mean_cov_high = mean(cov_high);
    result.mean_cov_total = mean(cov_total);
    result.mean_pareto_size = mean(pareto_sizes);
    result.archive_utilities = archive_utilities;

    fprintf('\n===== %s 结果 =====\n', variant_name);
    fprintf('  Utility: %.2f +/- %.2f\n', result.mean_utility, result.std_utility);
    fprintf('  Latency: %.2f s\n', result.mean_latency);
    fprintf('  Energy:  %.1f J\n', result.mean_energy);
    fprintf('  Pareto:  %.1f\n', result.mean_pareto_size);

    % 保存到临时文件
    save_file = fullfile(project_dir, 'experiments', sprintf('ablation_single_%s_%s.mat', variant_name, datestr(now, 'yyyymmdd_HHMMSS')));
    save(save_file, 'result');
    fprintf('结果已保存至: %s\n', save_file);
end

function pareto_front = archiveToFront(pareto_archive)
    pareto_front = zeros(length(pareto_archive), 3);
    for p_idx = 1:length(pareto_archive)
        pareto_front(p_idx, 1) = pareto_archive(p_idx).Utility;
        pareto_front(p_idx, 2) = pareto_archive(p_idx).Latency;
        pareto_front(p_idx, 3) = pareto_archive(p_idx).Energy;
    end
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
