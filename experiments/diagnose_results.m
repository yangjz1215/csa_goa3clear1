% diagnose_results.m - 从实验结果中提取关键指标进行诊断分析
clear; clc;

%% 1. 对比实验分析
comp_file = 'comparison_results_para_Map1_Medium_20260609_172450.mat';
data = load(comp_file);
results = data.results;

fprintf('========== 对比实验关键指标 ==========\n\n');

alg_names = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA2'};
for i = 1:length(alg_names)
    name = alg_names{i};
    if isfield(results, name)
        r = results.(name);
        fprintf('--- %s ---\n', name);
        fprintf('  Utility: %.2f ± %.2f\n', r.mean_utility, r.std_utility);
        fprintf('  Latency: %.4f ± %.4f\n', r.mean_latency, std(r.archive_latencies));
        fprintf('  Energy:  %.2f ± %.2f\n', r.mean_energy, std(r.archive_energies));
        fprintf('  HV(norm): %.4f ± %.4f\n', r.mean_hv_norm, r.std_hv_norm);
        fprintf('  IGD(norm): %.4f ± %.4f\n', r.mean_igd_norm, r.std_igd_norm);
        fprintf('  Spread: %.4f ± %.4f\n', r.mean_spread_norm, r.std_spread_norm);
        fprintf('  Runtime: %.2f s\n', r.mean_runtime);
        fprintf('  Pareto size: %.0f\n', r.mean_pareto_size);
        fprintf('  Success rate: %.4f\n', r.mean_success_rate);
        fprintf('\n');
    end
end

% 详细膝点分析
fprintf('\n========== 膝点详细对比 ==========\n');
for i = 1:length(alg_names)
    name = alg_names{i};
    if isfield(results, name) && isfield(results.(name), 'knee_utilities')
        r = results.(name);
        fprintf('%s: Knee U=%.2f L=%.4f E=%.2f\n', name, ...
            mean(r.knee_utilities), mean(r.knee_latencies), mean(r.knee_energies));
    end
end

% HV/IGD逐run对比
fprintf('\n========== HV逐run对比 ==========\n');
for i = 1:length(alg_names)
    name = alg_names{i};
    if isfield(results, name)
        hv = results.(name).hv_values;
        fprintf('%s: ', name);
        for j = 1:length(hv)
            if ~isnan(hv(j))
                fprintf('%.4f ', hv(j));
            end
        end
        fprintf('\n');
    end
end

fprintf('\n========== IGD逐run对比 ==========\n');
for i = 1:length(alg_names)
    name = alg_names{i};
    if isfield(results, name)
        igd = results.(name).igd_values;
        fprintf('%s: ', name);
        for j = 1:length(igd)
            if ~isnan(igd(j))
                fprintf('%.4f ', igd(j));
            end
        end
        fprintf('\n');
    end
end

% 综合评分（加入运行时间）
fprintf('\n========== 综合评分（HV+IGD+Runtime加权）==========\n');
hv_means = zeros(1, length(alg_names));
igd_means = zeros(1, length(alg_names));
runtime_means = zeros(1, length(alg_names));
for i = 1:length(alg_names)
    name = alg_names{i};
    if isfield(results, name)
        hv_means(i) = results.(name).mean_hv_norm;
        igd_means(i) = results.(name).mean_igd_norm;
        runtime_means(i) = results.(name).mean_runtime;
    end
end

% 归一化到[0,1]
hv_norm = (hv_means - min(hv_means)) / (max(hv_means) - min(hv_means) + 1e-9);
igd_norm = 1 - (igd_means - min(igd_means)) / (max(igd_means) - min(igd_means) + 1e-9); % IGD越小越好
time_norm = 1 - (runtime_means - min(runtime_means)) / (max(runtime_means) - min(runtime_means) + 1e-9); % 时间越短越好

for i = 1:length(alg_names)
    name = alg_names{i};
    score = 0.4 * hv_norm(i) + 0.4 * igd_norm(i) + 0.2 * time_norm(i);
    fprintf('%s: HV_score=%.3f IGD_score=%.3f Time_score=%.3f Combined=%.3f\n', ...
        name, hv_norm(i), igd_norm(i), time_norm(i), score);
end

%% 2. 消融实验分析
fprintf('\n\n========== 消融实验关键指标 ==========\n\n');
abl_file = 'ablation_results_para_Map1_Medium_20260609_202608.mat';
abl_data = load(abl_file);
abl_results = abl_data.results;

variant_names = {'proposed', 'no_subpop', 'no_goa_turn', 'no_goa_repulsion', 'no_pareto_leader', 'no_adaptive_weight'};
for i = 1:length(variant_names)
    name = variant_names{i};
    if isfield(abl_results, name)
        r = abl_results.(name);
        fprintf('--- %s ---\n', name);
        fprintf('  Utility: %.2f ± %.2f\n', r.mean_utility, r.std_utility);
        fprintf('  Latency: %.4f\n', r.mean_latency);
        fprintf('  Energy:  %.2f\n', r.mean_energy);
        fprintf('  HV(norm): %.4f ± %.4f\n', r.mean_hv_norm, r.std_hv_norm);
        fprintf('  IGD(norm): %.4f ± %.4f\n', r.mean_igd_norm, r.std_igd_norm);
        fprintf('  Spread: %.4f\n', r.mean_spread_norm);
        fprintf('  Pareto size: %.0f\n', r.mean_pareto_size);
        if isfield(r, 'knee_utilities')
            fprintf('  Knee: U=%.2f L=%.4f E=%.2f\n', r.mean_knee_utility, r.mean_knee_latency, r.mean_knee_energy);
        end
        fprintf('\n');
    end
end

% 消融对比：各变体相对于proposed的变化
fprintf('\n========== 消融对比（相对proposed变化%%）==========\n');
if isfield(abl_results, 'proposed')
    base = abl_results.proposed;
    for i = 2:length(variant_names)
        name = variant_names{i};
        if isfield(abl_results, name)
            r = abl_results.(name);
            util_change = (r.mean_utility - base.mean_utility) / base.mean_utility * 100;
            hv_change = (r.mean_hv_norm - base.mean_hv_norm) / abs(base.mean_hv_norm) * 100;
            igd_change = (r.mean_igd_norm - base.mean_igd_norm) / abs(base.mean_igd_norm) * 100;
            fprintf('%s: Utility%+.1f%% HV%+.1f%% IGD%+.1f%%\n', name, util_change, hv_change, igd_change);
        end
    end
end

fprintf('\n诊断完成。\n');
