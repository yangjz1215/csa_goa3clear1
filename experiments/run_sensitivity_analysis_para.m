function results = run_sensitivity_analysis_para(varargin)
    fprintf('========== 权重敏感性分析 (并行版) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'comparison_algorithms')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 10);
    addParameter(p, 'map_name', 'Map1_Medium');
    addParameter(p, 'n_workers', 4);
    parse(p, varargin{:});
    n_runs = p.Results.n_runs;
    map_name = p.Results.map_name;
    n_workers = p.Results.n_workers;

    maps_dir = fullfile(project_dir, 'maps');
    map_file = fullfile(maps_dir, [map_name, '.mat']);

    if ~exist(map_file, 'file')
        error('地图文件不存在: %s', map_file);
    end
    fprintf('--- 加载固定地图: %s ---\n', map_name);
    map_data = load(map_file);
    User = map_data.User;
    priorities = map_data.priorities;
    N_User = map_data.N_User;

    N_User_actual = size(User, 1);
    if N_User_actual == 200
        N_UAV = 8; config_suffix = 'Small';
    elseif N_User_actual == 500
        N_UAV = 15; config_suffix = 'Medium';
    else
        N_UAV = 25; config_suffix = 'Large';
    end
    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);

    config_file = fullfile(maps_dir, ['Map_', config_suffix, '_Config.mat']);
    if exist(config_file, 'file')
        config_data = load(config_file);
        RRH = config_data.RRH;
        RRH_type = config_data.RRH_type;
        N_RRH = config_data.N_RRH;
        N_eRRH = config_data.N_eRRH;
        UAV_type = config_data.UAV_type;
        params.D = config_data.D;
        params.C = config_data.C;
        params.DT = config_data.DT;
    else
        Lb = [0, 0]; Ub = [1000, 1000];
        N_RRH = 10; N_eRRH = 4;
        RRH = Lb + (Ub - Lb) .* rand(N_RRH, 2);
        RRH_type = zeros(N_RRH, 1); RRH_type(1:N_eRRH) = 1;
        UAV_type = zeros(N_UAV, 1); UAV_type(1:floor(N_UAV*0.3)) = 1;
        params.D = []; params.C = []; params.DT = [];
    end

    Lb = [0, 0]; Ub = [1000, 1000];

    base_params = struct();
    base_params.RRH = RRH;

    if exist('config_data', 'var') && isfield(config_data, 'D') && ~isempty(config_data.D)
        base_params.D = config_data.D;
        base_params.C = config_data.C;
    else
        base_params.D = ones(N_User, 1) * 2e6;
        base_params.C = ones(N_User, 1) * 0.5e9;
    end

    base_params.D_UU = 10; base_params.D_RU = 10;
    base_params.cover_radius = 150; base_params.RRH_radius = 150;
    base_params.E_max = 50000; base_params.k_move = 15;
    base_params.Pho = 100; base_params.ki = 1e-27;
    base_params.PtxU = 0.5; base_params.PtxEU = 5;
    base_params.Ptx = 1; base_params.PtxR = 10;
    base_params.alpha0 = 1.42e-4; base_params.sigma2 = 3.98e-12;
    base_params.f_eUAV = 2e9; base_params.f_eRRH = 2e9; base_params.f_BBU = 50e9;
    base_params.F_total = 10e9; base_params.B_total = 10e6;
    base_params.B_total_relay = 20e6;
    base_params.max_latency = 1.0; base_params.noise = 1e-13;
    base_params.kappa = 1e-28; base_params.P_tx = 0.5;
    base_params.FES_max = 150; base_params.K = 40;
    base_params.enable_early_stop = false; base_params.enable_smart_stop = false;
    base_params.enable_bilevel = true;
    base_params.enable_pareto_leader = true;
    base_params.G_weights = [0.4, 0.3, 0.3];
    base_params.subpop_params = struct();
    base_params.subpop_params.mu0 = [500, 500];
    base_params.subpop_params.sigma0 = [150, 150; 120, 120; 80, 80];
    base_params.subpop_params.sigma_min = [5, 8, 3];
    base_params.subpop_params.q = [0.6, 0.5, 0.4];
    base_params.subpop_params.beta = [0.8, 0.7, 0.6];

    %% ========================================================================
    %% PHASE 1: 算法底层超参数 (K, q) 网格搜索与调优
    %% ========================================================================
    fprintf('\n>>> [Phase 1] 开始底层超参数网格搜索...\n');

    K_list = [10, 15, 20, 40];
    q_list = [0.4, 0.6, 0.8];
    [K_grid, q_grid] = meshgrid(K_list, q_list);
    num_algo_configs = numel(K_grid);
    algo_hv_results = zeros(num_algo_configs, 1);

    parfor ac_idx = 1:num_algo_configs
        test_K = K_grid(ac_idx);
        test_q = q_grid(ac_idx);
        fprintf('  [Phase 1] K=%d, q=%.1f\n', test_K, test_q);

        p_struct = base_params;
        p_struct.K = test_K;
        p_struct.subpop_params.q = [test_q, test_q, test_q];

        [~, ~, ~, ~, ~, ~, ~, ~, ~, pareto_archive] = ...
            cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, ...
            Ub, Lb, p_struct, priorities);

        algo_hv_results(ac_idx) = 0;
        if ~isempty(pareto_archive)
            arch_U = [pareto_archive.Utility];
            arch_L = [pareto_archive.Latency];
            arch_E = [pareto_archive.Energy];
            max_U = sum(priorities);
            max_L = N_User * base_params.max_latency;
            max_E = base_params.E_max * N_UAV;
            n_U = (max_U - arch_U) / (max_U + 1e-6);
            n_L = arch_L / (max_L + 1e-6);
            n_E = arch_E / (max_E + 1e-6);
            normalized_objs = [n_U', n_L', n_E'];
            algo_hv_results(ac_idx) = hypervolume(normalized_objs, [1.1, 1.1, 1.1]);
        end
    end

    [best_hv_algo, best_algo_idx] = max(algo_hv_results);
    best_K = K_grid(best_algo_idx);
    best_q = q_grid(best_algo_idx);

    fprintf('\n>>> [Phase 1 完成] 最优超参数: K=%d, q=%.1f (HV=%.4f)\n', best_K, best_q, best_hv_algo);

    base_params.K = best_K;
    base_params.subpop_params.q = [best_q, best_q, best_q];
    fprintf('    [自动注入] 已将最优 K/q 写入 base_params，后续 Phase 2 将使用此配置\n');

    save(fullfile(project_dir, 'experiments', 'best_algo_hyperparams.mat'), 'best_K', 'best_q');

    %% ========================================================================
    %% PHASE 2: 决策权重偏好 (w) 敏感性分析 (使用最优 K 和 q)
    %% ========================================================================
    fprintf('\n>>> [Phase 2] 开始权重偏好敏感性分析...\n');

    W0 = [0.60, 0.10, 0.30;
          0.40, 0.40, 0.20;
          0.30, 0.20, 0.50];

    weight_configs = cell(10, 1);
    weight_configs{1} = W0;
    weight_configs{2} = [0.66, 0.07, 0.27; W0(2,:); W0(3,:)];
    weight_configs{3} = [0.54, 0.13, 0.33; W0(2,:); W0(3,:)];
    weight_configs{4} = [W0(1,:); 0.36, 0.44, 0.20; W0(3,:)];
    weight_configs{5} = [W0(1,:); 0.44, 0.36, 0.20; W0(3,:)];
    weight_configs{6} = [W0(1,:); W0(2,:); 0.27, 0.18, 0.55];
    weight_configs{7} = [W0(1,:); W0(2,:); 0.33, 0.22, 0.45];
    weight_configs{8} = [0.65, 0.10, 0.25; 0.35, 0.45, 0.20; 0.25, 0.20, 0.55];
    weight_configs{9} = [0.55, 0.15, 0.30; 0.45, 0.35, 0.20; 0.35, 0.25, 0.40];
    weight_configs{10} = [0.33, 0.33, 0.34; 0.33, 0.34, 0.33; 0.34, 0.33, 0.33];

    config_names = {
        '1. Baseline (基准)', ...
        '2. G1 Utility +10%', ...
        '3. G1 Utility -10%', ...
        '4. G2 Latency +10%', ...
        '5. G2 Latency -10%', ...
        '6. G3 Energy +10%', ...
        '7. G3 Energy -10%', ...
        '8. All Stronger', ...
        '9. All Weaker', ...
        '10. Uniform (对照)'
    };

    num_configs = length(weight_configs);

    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local', n_workers);
    end

    hv_results = zeros(num_configs, n_runs);
    archive_sizes = zeros(num_configs, n_runs);
    best_utils = zeros(num_configs, n_runs);
    best_lats = zeros(num_configs, n_runs);
    best_nrgs = zeros(num_configs, n_runs);

    fprintf('开始并行运算，共 %d 组实验，每组 %d 次运行...\n', num_configs, n_runs);

    parfor cfg_idx = 1:num_configs
        cfg_hvs = zeros(1, n_runs);
        cfg_sizes = zeros(1, n_runs);
        cfg_utils = zeros(1, n_runs);
        cfg_lats = zeros(1, n_runs);
        cfg_nrgs = zeros(1, n_runs);

        test_weights = weight_configs{cfg_idx};

        for r = 1:n_runs
            params_r = base_params;
            params_r.test_weights = test_weights;

            [~, ~, ~, ~, ~, ~, ~, ~, ~, pareto_archive] = ...
                cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, ...
                Ub, Lb, params_r, priorities);

            if ~isempty(pareto_archive)
                cfg_sizes(r) = length(pareto_archive);
                arch_util = [pareto_archive.Utility];
                arch_lat = [pareto_archive.Latency];
                arch_energy = [pareto_archive.Energy];

                cfg_utils(r) = max(arch_util);
                cfg_lats(r) = min(arch_lat);
                cfg_nrgs(r) = min(arch_energy);

                max_U = sum(priorities);
                max_L = N_User * base_params.max_latency;
                max_E = base_params.E_max * N_UAV;

                n_U = (max_U - arch_util) / (max_U + 1e-6);
                n_L = arch_lat / (max_L + 1e-6);
                n_E = arch_energy / (max_E + 1e-6);
                normalized_objs = [n_U', n_L', n_E'];

                cfg_hvs(r) = hypervolume(normalized_objs, [1.1, 1.1, 1.1]);
            else
                cfg_sizes(r) = 0; cfg_hvs(r) = 0;
                cfg_utils(r) = 0; cfg_lats(r) = Inf; cfg_nrgs(r) = Inf;
            end
        end

        hv_results(cfg_idx, :) = cfg_hvs;
        archive_sizes(cfg_idx, :) = cfg_sizes;
        best_utils(cfg_idx, :) = cfg_utils;
        best_lats(cfg_idx, :) = cfg_lats;
        best_nrgs(cfg_idx, :) = cfg_nrgs;

        fprintf('  配置 %d/%d (%s) 完成\n', cfg_idx, num_configs, config_names{cfg_idx});
    end

    mean_hv = mean(hv_results, 2);
    std_hv = std(hv_results, 0, 2);
    mean_size = mean(archive_sizes, 2);
    mean_util = mean(best_utils, 2);
    mean_lat = mean(best_lats, 2);
    mean_nrg = mean(best_nrgs, 2);

    fprintf('\n========== 权重敏感性分析结果 ==========\n');
    fprintf('%-30s %8s %10s %10s %10s %10s\n', '配置', 'HV均值', 'HV标准差', '最大效用', '最小能耗', '解数量');
    fprintf('--------------------------------------------------------------------------------\n');
    for i = 1:num_configs
        fprintf('%-30s %8.4f %10.4f %10.1f %10.1f %10.1f\n', ...
            config_names{i}, mean_hv(i), std_hv(i), mean_util(i), mean_nrg(i), mean_size(i));
    end

    T = table(config_names', mean_hv, std_hv, mean_util, mean_lat, mean_nrg, mean_size, ...
        'VariableNames', {'Configuration', 'HV_Mean', 'HV_Std', 'Max_Utility', 'Min_Latency', 'Min_Energy', 'ArchiveSize'});

    csv_file = fullfile(project_dir, 'experiments', sprintf('weight_sensitivity_%s.csv', datestr(now, 'yyyymmdd_HHMMSS')));
    writetable(T, csv_file);
    fprintf('\n结果已保存至: %s\n', csv_file);

    baseline_hv = mean_hv(1);
    fprintf('\n========== 相对基准波动分析 ==========\n');
    fprintf('%-30s %10s\n', '配置', '相对基准HV偏差');
    fprintf('--------------------------------\n');
    for i = 1:num_configs
        pct = (mean_hv(i) - baseline_hv) / baseline_hv * 100;
        fprintf('%-30s %8.2f%%\n', config_names{i}, pct);
    end

    %% ========================================================================
    %% PHASE 3: Tier 2 核心超参数 (phi_t 权重 + beta) 敏感性分析
    %% ========================================================================
    fprintf('\n>>> [Phase 3] 开始 Tier 2 核心超参数敏感性分析 (±10%%)...\n');

    tier2_num_configs = 12;
    tier2_config_names = cell(tier2_num_configs, 1);
    tier2_inject_field = cell(tier2_num_configs, 1);
    tier2_inject_value = zeros(tier2_num_configs, 1);

    % phi_t w1 (progress): [0.252, 0.280, 0.308]
    tier2_config_names{1}  = 'phi_t w1 (progress) -10%';  tier2_inject_field{1} = 'phase_w_progress'; tier2_inject_value(1) = 0.252;
    tier2_config_names{2}  = 'phi_t w1 (progress) 基准';   tier2_inject_field{2} = 'phase_w_progress'; tier2_inject_value(2) = 0.280;
    tier2_config_names{3}  = 'phi_t w1 (progress) +10%';  tier2_inject_field{3} = 'phase_w_progress'; tier2_inject_value(3) = 0.308;
    % phi_t w2 (coverage): [0.378, 0.420, 0.462]
    tier2_config_names{4}  = 'phi_t w2 (coverage) -10%';  tier2_inject_field{4} = 'phase_w_cov';      tier2_inject_value(4) = 0.378;
    tier2_config_names{5}  = 'phi_t w2 (coverage) 基准';   tier2_inject_field{5} = 'phase_w_cov';      tier2_inject_value(5) = 0.420;
    tier2_config_names{6}  = 'phi_t w2 (coverage) +10%';  tier2_inject_field{6} = 'phase_w_cov';      tier2_inject_value(6) = 0.462;
    % phi_t w3 (inner): [0.270, 0.300, 0.330]
    tier2_config_names{7}  = 'phi_t w3 (inner) -10%';     tier2_inject_field{7} = 'phase_w_inner';    tier2_inject_value(7) = 0.270;
    tier2_config_names{8}  = 'phi_t w3 (inner) 基准';      tier2_inject_field{8} = 'phase_w_inner';    tier2_inject_value(8) = 0.300;
    tier2_config_names{9}  = 'phi_t w3 (inner) +10%';     tier2_inject_field{9} = 'phase_w_inner';    tier2_inject_value(9) = 0.330;
    % beta (uniform): [0.63, 0.70, 0.77]
    tier2_config_names{10} = 'beta (uniform) -10%';       tier2_inject_field{10} = 'beta';            tier2_inject_value(10) = 0.63;
    tier2_config_names{11} = 'beta (uniform) 基准';        tier2_inject_field{11} = 'beta';            tier2_inject_value(11) = 0.70;
    tier2_config_names{12} = 'beta (uniform) +10%';       tier2_inject_field{12} = 'beta';            tier2_inject_value(12) = 0.77;

    tier2_hv_results = zeros(tier2_num_configs, n_runs);
    tier2_archive_sizes = zeros(tier2_num_configs, n_runs);
    tier2_best_utils = zeros(tier2_num_configs, n_runs);
    tier2_best_lats = zeros(tier2_num_configs, n_runs);
    tier2_best_nrgs = zeros(tier2_num_configs, n_runs);

    fprintf('开始并行运算，共 %d 组 Tier 2 实验，每组 %d 次运行...\n', tier2_num_configs, n_runs);

    parfor cfg_idx = 1:tier2_num_configs
        cfg_hvs = zeros(1, n_runs);
        cfg_sizes = zeros(1, n_runs);
        cfg_utils = zeros(1, n_runs);
        cfg_lats = zeros(1, n_runs);
        cfg_nrgs = zeros(1, n_runs);

        for r = 1:n_runs
            params_r = base_params;
            field_name = tier2_inject_field{cfg_idx};
            val = tier2_inject_value(cfg_idx);

            % phi_t 权重通过 params 字段注入，computePhasePhi.m 自动归一化
            if strcmp(field_name, 'phase_w_progress') || strcmp(field_name, 'phase_w_cov') || strcmp(field_name, 'phase_w_inner')
                params_r.(field_name) = val;
            elseif strcmp(field_name, 'beta')
                params_r.subpop_params.beta = [val, val, val];
            end

            [~, ~, ~, ~, ~, ~, ~, ~, ~, pareto_archive] = ...
                cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, ...
                Ub, Lb, params_r, priorities);

            if ~isempty(pareto_archive)
                cfg_sizes(r) = length(pareto_archive);
                arch_util = [pareto_archive.Utility];
                arch_lat = [pareto_archive.Latency];
                arch_energy = [pareto_archive.Energy];

                cfg_utils(r) = max(arch_util);
                cfg_lats(r) = min(arch_lat);
                cfg_nrgs(r) = min(arch_energy);

                max_U = sum(priorities);
                max_L = N_User * base_params.max_latency;
                max_E = base_params.E_max * N_UAV;

                n_U = (max_U - arch_util) / (max_U + 1e-6);
                n_L = arch_lat / (max_L + 1e-6);
                n_E = arch_energy / (max_E + 1e-6);
                normalized_objs = [n_U', n_L', n_E'];

                cfg_hvs(r) = hypervolume(normalized_objs, [1.1, 1.1, 1.1]);
            else
                cfg_sizes(r) = 0; cfg_hvs(r) = 0;
                cfg_utils(r) = 0; cfg_lats(r) = Inf; cfg_nrgs(r) = Inf;
            end
        end

        tier2_hv_results(cfg_idx, :) = cfg_hvs;
        tier2_archive_sizes(cfg_idx, :) = cfg_sizes;
        tier2_best_utils(cfg_idx, :) = cfg_utils;
        tier2_best_lats(cfg_idx, :) = cfg_lats;
        tier2_best_nrgs(cfg_idx, :) = cfg_nrgs;

        fprintf('  Tier2 配置 %d/%d (%s) 完成\n', cfg_idx, tier2_num_configs, tier2_config_names{cfg_idx});
    end

    tier2_mean_hv = mean(tier2_hv_results, 2);
    tier2_std_hv = std(tier2_hv_results, 0, 2);
    tier2_mean_size = mean(tier2_archive_sizes, 2);
    tier2_mean_util = mean(tier2_best_utils, 2);
    tier2_mean_lat = mean(tier2_best_lats, 2);
    tier2_mean_nrg = mean(tier2_best_nrgs, 2);

    fprintf('\n========== Tier 2 超参数敏感性分析结果 ==========\n');
    fprintf('%-35s %8s %10s %10s %10s %10s\n', '配置', 'HV均值', 'HV标准差', '最大效用', '最小能耗', '解数量');
    fprintf('------------------------------------------------------------------------------------\n');
    for i = 1:tier2_num_configs
        fprintf('%-35s %8.4f %10.4f %10.1f %10.1f %10.1f\n', ...
            tier2_config_names{i}, tier2_mean_hv(i), tier2_std_hv(i), tier2_mean_util(i), tier2_mean_nrg(i), tier2_mean_size(i));
    end

    % 自动选出 Tier 2 最优参数并注入 base_params
    [best_tier2_hv, best_tier2_idx] = max(tier2_mean_hv);
    best_phase_w_progress = 0.28;
    best_phase_w_cov = 0.42;
    best_phase_w_inner = 0.30;
    best_beta = 0.70;

    field_name_best = tier2_inject_field{best_tier2_idx};
    val_best = tier2_inject_value(best_tier2_idx);
    if strcmp(field_name_best, 'phase_w_progress')
        best_phase_w_progress = val_best;
    elseif strcmp(field_name_best, 'phase_w_cov')
        best_phase_w_cov = val_best;
    elseif strcmp(field_name_best, 'phase_w_inner')
        best_phase_w_inner = val_best;
    elseif strcmp(field_name_best, 'beta')
        best_beta = val_best;
    end

    base_params.phase_w_progress = best_phase_w_progress;
    base_params.phase_w_cov = best_phase_w_cov;
    base_params.phase_w_inner = best_phase_w_inner;
    base_params.subpop_params.beta = [best_beta, best_beta, best_beta];

    fprintf('\n>>> [Phase 3 完成] 最优 Tier 2 超参数:\n');
    fprintf('    phi_t 权重: progress=%.3f, coverage=%.3f, inner=%.3f\n', best_phase_w_progress, best_phase_w_cov, best_phase_w_inner);
    fprintf('    beta=%.2f (HV=%.4f)\n', best_beta, best_tier2_hv);
    fprintf('    [自动注入] 已将最优 phi_t 权重和 beta 写入 base_params\n');

    %% 保存所有结果
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_filename = sprintf('sensitivity_results_para_%s_%s.mat', map_name, timestamp);
    save_path = fullfile(project_dir, 'experiments', save_filename);
    params = base_params;
    save(save_path, ...
        'config_names', 'weight_configs', 'archive_sizes', 'best_utils', 'hv_results', ...
        'map_name', 'params', 'mean_hv', 'std_hv', 'mean_util', 'mean_lat', 'mean_nrg', 'mean_size', ...
        'K_grid', 'q_grid', 'algo_hv_results', 'best_K', 'best_q', 'best_hv_algo', ...
        'tier2_config_names', 'tier2_hv_results', 'tier2_mean_hv', 'tier2_std_hv', ...
        'tier2_mean_util', 'tier2_mean_lat', 'tier2_mean_nrg', 'tier2_mean_size', ...
        'tier2_inject_field', 'tier2_inject_value', ...
        'best_phase_w_progress', 'best_phase_w_cov', 'best_phase_w_inner', 'best_beta');
    fprintf('敏感性分析核心数据已成功保存至: %s\n', save_path);

    % 保存最优超参数到 best_algo_hyperparams.mat（追加新变量，保留已有 best_K/best_q）
    hyperparam_save = fullfile(project_dir, 'experiments', 'best_algo_hyperparams.mat');
    if exist(hyperparam_save, 'file')
        existing = load(hyperparam_save);
        best_K = existing.best_K;
        best_q = existing.best_q;
    end
    save(hyperparam_save, 'best_K', 'best_q', 'best_phase_w_progress', 'best_phase_w_cov', 'best_phase_w_inner', 'best_beta');
    fprintf('最优超参数已保存至: %s\n', hyperparam_save);
    fprintf('>> 提示: 请前往 visualization 文件夹运行 plot_sensitivity.m 生成最终的图表和表格\n');

    results = struct();
    results.hv_mean = mean_hv; results.hv_std = std_hv;
    results.config_names = config_names;
    results.table = T;
    results.tier2_hv_mean = tier2_mean_hv; results.tier2_hv_std = tier2_std_hv;
    results.tier2_config_names = tier2_config_names;
    fprintf('\n========== 超参数敏感性分析全部完成 ==========\n');
end