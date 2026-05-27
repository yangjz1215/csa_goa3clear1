function plot_comparison_all(data_file, output_dir)
% plot_comparison_all - 对比实验完整可视化（收敛曲线/Pareto散点/箱型图）
% 支持自动识别场景（Small/Medium/Large）并输出到对应子文件夹
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'comparison_results_para_*.mat');
    files = dir(data_file);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        data_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No comparison results file found');
    end
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures', 'comparison');
end

[~, base_name, ~] = fileparts(data_file);
scene = 'Unknown';
if contains(base_name, 'Small', 'IgnoreCase', true)
    scene = 'Small';
elseif contains(base_name, 'Medium', 'IgnoreCase', true)
    scene = 'Medium';
elseif contains(base_name, 'Large', 'IgnoreCase', true)
    scene = 'Large';
end

output_dir = fullfile(output_dir, scene);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

fig1 = figure('Position', [100, 100, 900, 500]);
set(gcf, 'Color', 'w');

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

algorithms = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA2'};
labels = {'cSA-GOA (Proposed)', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA-II'};
marker_styles = {'*', 'o', 's', '^', 'd', 'v'};
sizes = [80, 40, 40, 40, 40, 40];
n_algs = length(algorithms);

subplot(1, 2, 1);
hold on;
target_len = 0;
if exist('params', 'var') && isfield(params, 'FES_max') && ~isempty(params.FES_max)
    target_len = params.FES_max;
end
if target_len <= 0
    for a = 1:length(algorithms)
        alg = algorithms{a};
        if isfield(results, alg) && isfield(results.(alg), 'convergence_curves')
            curves = results.(alg).convergence_curves;
            if iscell(curves) && ~isempty(curves)
                target_len = max(target_len, max(cellfun(@length, curves)));
            end
        end
    end
end
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'convergence_curves')
        curves = results.(alg).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            curve_matrix = nan(min(length(curves), 30), target_len);
            row_count = 0;
            for run = 1:min(length(curves), 30)
                if ~isempty(curves{run})
                    c = double(curves{run}(:)');
                    c = cummax(c);
                    if length(c) >= target_len
                        c_plot = c(1:target_len);
                    else
                        c_plot = [c, repmat(c(end), 1, target_len - length(c))];
                    end
                    row_count = row_count + 1;
                    curve_matrix(row_count, :) = c_plot;
                end
            end
            if row_count > 0
                avg_curve = mean(curve_matrix(1:row_count, :), 1, 'omitnan');
                if a == 1
                    lw = 2.5;
                else
                    lw = 1.8;
                end
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors(a,:), 'LineWidth', lw, 'DisplayName', labels{a});
            end
        end
    end
end
xlabel('Generations', 'FontWeight', 'bold');
ylabel('Best-so-far Utility (Priority Sum) \uparrow', 'FontWeight', 'bold');
title(['Utility Convergence (', scene, ')'], 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
xlim([1, target_len]);
if exist('map_data', 'var') && isfield(map_data, 'priorities') && ~isempty(map_data.priorities)
    ylim([0, sum(map_data.priorities) * 1.02]);
end
grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

subplot(1, 2, 2);
hold on;
hv_means = zeros(1, length(algorithms));
hv_stds = zeros(1, length(algorithms));
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_means(a) = mean(vals);
        hv_stds(a) = std(vals);
    end
end
errorbar(1:length(algorithms), hv_means, hv_stds, '-o', 'Color', [0.0000, 0.4470, 0.7410], ...
    'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.0000, 0.4470, 0.7410], 'Capsize', 6);
xlabel('Algorithm', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title(['HV Comparison (Mean ± Std) - ', scene], 'FontWeight', 'bold');
xticks(1:length(algorithms));
xticklabels(labels);
xtickangle(45);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

saveas(fig1, fullfile(output_dir, 'comparison_convergence.fig'));
saveas(fig1, fullfile(output_dir, 'comparison_convergence.png'));
close(fig1);

fig2 = figure('Position', [100, 100, 800, 600]);
set(gcf, 'Color', 'w');
hold on;
grid on;
view(45, 30);

sz = 60;
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'pareto_fronts')
        pfs = results.(alg).pareto_fronts;

        best_run = 1; max_size = 0;
        for run = 1:length(pfs)
            if ~isempty(pfs{run}) && size(pfs{run},1) > max_size
                max_size = size(pfs{run},1);
                best_run = run;
            end
        end

        best_pf = pfs{best_run};
        if ~isempty(best_pf) && size(best_pf,2) == 3
            scatter3(best_pf(:,1), best_pf(:,2), best_pf(:,3), sz, ...
                'MarkerEdgeColor', colors(a,:), 'MarkerFaceColor', colors(a,:), ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', labels{a});
        end
    end
end

xlabel('System Utility (Priority Sum)', 'FontWeight', 'bold');
ylabel('Total Latency (s)', 'FontWeight', 'bold');
zlabel('Total Energy (J)', 'FontWeight', 'bold');
title(['3D Pareto Front Comparison (', scene, ')'], 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 10);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

saveas(fig2, fullfile(output_dir, 'comparison_pareto_3D.fig'));
saveas(fig2, fullfile(output_dir, 'comparison_pareto_3D.png'));
close(fig2);

fig3 = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

hv_all = [];
group_idx = [];
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_all = [hv_all; vals(:)];
        group_idx = [group_idx; a * ones(length(vals(:)), 1)];
    end
end

box_colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

hold on;
for a = 1:length(algorithms)
    idx = (group_idx == a);
    if any(idx)
        b = boxchart(group_idx(idx), hv_all(idx));
        b.BoxFaceColor = box_colors(a, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:length(algorithms));
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title(['HV Distribution (30 Independent Runs) - ', scene], 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3;
box on;

saveas(fig3, fullfile(output_dir, 'comparison_boxplot.fig'));
saveas(fig3, fullfile(output_dir, 'comparison_boxplot.png'));
close(fig3);

fprintf('Comparison charts saved to %s\n', output_dir);
fprintf('  - comparison_convergence.fig/png (收敛曲线 + HV误差棒)\n');
fprintf('  - comparison_pareto_3D.fig/png (3D Pareto 前沿)\n');
fprintf('  - comparison_boxplot.fig/png (HV分布箱型图)\n');

fig4 = figure('Name', 'Pareto Projections', 'Position', [100, 100, 1500, 450]);
set(gcf, 'Color', 'w');

proj_configs = {
    2, 1, 'Total Latency (s) \downarrow', 'System Utility \uparrow', 1;
    3, 1, 'Total Energy (J) \downarrow', 'System Utility \uparrow', 2;
    2, 3, 'Total Latency (s) \downarrow', 'Total Energy (J) \downarrow', 3
};

for p = 1:3
    subplot(1, 3, proj_configs{p, 5});
    hold on; grid on;

    for a = 1:length(algorithms)
        alg = algorithms{a};
        if ~isfield(results, alg) || ~isfield(results.(alg), 'pareto_fronts'), continue; end

        pfs = results.(alg).pareto_fronts;
        all_data = [];
        for run = 1:length(pfs)
            if ~isempty(pfs{run}), all_data = [all_data; pfs{run}]; end
        end

        if ~isempty(all_data)
            if strcmp(alg, 'cSA_GOA')
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    45, colors(a,:), 'o', 'filled', ...
                    'MarkerFaceAlpha', 1.0, ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], ...
                    'LineWidth', 0.5, ...
                    'DisplayName', labels{a});
            else
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    20, colors(a,:), 'o', 'filled', ...
                    'MarkerFaceAlpha', 0.08, ...
                    'MarkerEdgeColor', 'none', ...
                    'HandleVisibility', 'on', ...
                    'DisplayName', labels{a});
            end
        end
    end

    xlabel(proj_configs{p, 3}, 'FontWeight', 'bold');
    ylabel(proj_configs{p, 4}, 'FontWeight', 'bold');
    title(['Projection: ', proj_configs{p, 3}, ' vs ', proj_configs{p, 4}]);

    if p == 1, legend('Location', 'best', 'FontSize', 7); end
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.1);
    box on;
end

saveas(fig4, fullfile(output_dir, 'comparison_projections_2d.png'));
saveas(fig4, fullfile(output_dir, 'comparison_projections_2d.fig'));
close(fig4);

num_algs = length(algorithms);

% 检查旧 mat 是否有 runtime/success_rate 数据
has_runtimes = false;
for a = 1:num_algs
    if isfield(results, algorithms{a}) && isfield(results.(algorithms{a}), 'runtimes')
        has_runtimes = true; break;
    end
end

if has_runtimes
fprintf('正在绘制收敛速度与成功率图表...\n');

n_runs_actual = 0;
for a = 1:num_algs
    if isfield(results, algorithms{a}) && isfield(results.(algorithms{a}), 'runtimes')
        n_runs_actual = max(n_runs_actual, length(results.(algorithms{a}).runtimes));
    end
end
all_runtimes = nan(n_runs_actual, num_algs);
all_success = nan(n_runs_actual, num_algs);

for a = 1:num_algs
    alg_name = algorithms{a};
    if isfield(results, alg_name) && isfield(results.(alg_name), 'runtimes')
        rt = results.(alg_name).runtimes(:);
        all_runtimes(1:length(rt), a) = rt;
        if isfield(results.(alg_name), 'success_rates')
            sr = results.(alg_name).success_rates(:) * 100;
            all_success(1:length(sr), a) = sr;
        end
    end
end

fig_runtime = figure('Position', [100, 100, 600, 450]);
set(gcf, 'Color', 'w');
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1.2);

total_iterations = 300;
if exist('params', 'var') && isfield(params, 'FES_max') && ~isempty(params.FES_max)
    total_iterations = params.FES_max;
end
mean_rt = mean(all_runtimes, 1, 'omitnan') / total_iterations;
std_rt = std(all_runtimes, 0, 1, 'omitnan') / total_iterations;

for a = 1:num_algs
    if ~isnan(mean_rt(a))
        bar(a, mean_rt(a), 0.6, 'FaceColor', colors(a,:), 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
    end
end
errorbar(1:num_algs, mean_rt, std_rt, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
ylabel('Average Runtime per Iteration (Seconds)', 'FontWeight', 'bold', 'FontSize', 12);
title(['Algorithm Runtime Comparison (', scene, ')'], 'FontWeight', 'bold', 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

exportgraphics(fig_runtime, fullfile(output_dir, 'comparison_runtime_bar.png'), 'Resolution', 300);
saveas(fig_runtime, fullfile(output_dir, 'comparison_runtime_bar.fig'));
close(fig_runtime);

fig_success = figure('Position', [150, 150, 600, 450]);
set(gcf, 'Color', 'w');

hold on;
for a = 1:num_algs
    valid = ~isnan(all_success(:, a));
    if any(valid)
        b = boxchart(a * ones(sum(valid), 1), all_success(valid, a));
        b.BoxFaceColor = colors(a, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end
set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
ylabel('Task Execution Success Rate (%)', 'FontWeight', 'bold', 'FontSize', 12);
title(['Task Success Rate Distribution (30 Independent Runs) - ', scene], 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

exportgraphics(fig_success, fullfile(output_dir, 'comparison_success_rate_boxplot.png'), 'Resolution', 300);
saveas(fig_success, fullfile(output_dir, 'comparison_success_rate_boxplot.fig'));
close(fig_success);

fprintf('Runtime and success rate charts saved to %s\n', output_dir);
else
    fprintf('Old mat file detected - skipping runtime/success_rate charts (no data)\n');
end

% 检查旧 mat 是否有 IGD/GD/Spread 数据
has_igd = false;
for a = 1:num_algs
    if isfield(results, algorithms{a}) && isfield(results.(algorithms{a}), 'igd_values')
        has_igd = true; break;
    end
end

if has_igd
fprintf('正在绘制 IGD、GD 与 Spread 箱线图...\n');

n_runs_metrics = 0;
for a = 1:num_algs
    if isfield(results, algorithms{a}) && isfield(results.(algorithms{a}), 'igd_values')
        n_runs_metrics = max(n_runs_metrics, length(results.(algorithms{a}).igd_values));
    end
end
all_igd = nan(n_runs_metrics, num_algs);
all_gd = nan(n_runs_metrics, num_algs);
all_spread = nan(n_runs_metrics, num_algs);

for a = 1:num_algs
    alg_name = algorithms{a};
    if ~isfield(results, alg_name), continue; end
    if isfield(results.(alg_name), 'igd_values')
        temp_igd = results.(alg_name).igd_values(:);
        all_igd(1:length(temp_igd), a) = temp_igd;
    end
    if isfield(results.(alg_name), 'gd_values')
        gd = results.(alg_name).gd_values(:);
        all_gd(1:length(gd), a) = gd;
    end
    if isfield(results.(alg_name), 'spread_values')
        sp = results.(alg_name).spread_values(:);
        all_spread(1:length(sp), a) = sp;
    end
end

fig_igd = figure('Position', [200, 200, 600, 450]);
set(gcf, 'Color', 'w');
hold on;
for a = 1:num_algs
    b = boxchart(a * ones(size(all_igd, 1), 1), all_igd(:, a));
    b.BoxFaceColor = colors(a, :);
    b.BoxFaceAlpha = 0.6;
    b.MarkerStyle = 'o';
    b.MarkerColor = [0.2, 0.2, 0.2];
    b.LineWidth = 1.5;
end
set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
ylabel('Inverted Generational Distance (IGD) \downarrow', 'FontWeight', 'bold', 'FontSize', 12);
title(['Convergence Evaluation: IGD Metric (Lower is Better) - ', scene], 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;
exportgraphics(fig_igd, fullfile(output_dir, 'comparison_igd_boxplot.png'), 'Resolution', 300);
saveas(fig_igd, fullfile(output_dir, 'comparison_igd_boxplot.fig'));
close(fig_igd);

if any(all_gd(~isnan(all_gd)) > 0)
    fig_gd = figure('Position', [220, 220, 600, 450]);
    set(gcf, 'Color', 'w');
    hold on;
    for a = 1:num_algs
        b = boxchart(a * ones(size(all_gd, 1), 1), all_gd(:, a));
        b.BoxFaceColor = colors(a, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
    set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
    ylabel('Generational Distance (GD) \downarrow', 'FontWeight', 'bold', 'FontSize', 12);
    title(['Convergence Evaluation: GD Metric (Lower is Better) - ', scene], 'FontWeight', 'bold', 'FontSize', 14);
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1.2);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
    box on;
    exportgraphics(fig_gd, fullfile(output_dir, 'comparison_gd_boxplot.png'), 'Resolution', 300);
    saveas(fig_gd, fullfile(output_dir, 'comparison_gd_boxplot.fig'));
    close(fig_gd);
end

fig_spread = figure('Position', [250, 250, 600, 450]);
set(gcf, 'Color', 'w');
hold on;
for a = 1:num_algs
    b = boxchart(a * ones(size(all_spread, 1), 1), all_spread(:, a));
    b.BoxFaceColor = colors(a, :);
    b.BoxFaceAlpha = 0.6;
    b.MarkerStyle = 'o';
    b.MarkerColor = [0.2, 0.2, 0.2];
    b.LineWidth = 1.5;
end
set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
ylabel('Spread Metric \uparrow', 'FontWeight', 'bold', 'FontSize', 12);
title(['Diversity Evaluation: Spread Metric (Higher is Better) - ', scene], 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;
exportgraphics(fig_spread, fullfile(output_dir, 'comparison_spread_boxplot.png'), 'Resolution', 300);
saveas(fig_spread, fullfile(output_dir, 'comparison_spread_boxplot.fig'));
close(fig_spread);

fprintf('IGD and Spread charts saved to %s\n', output_dir);
else
    fprintf('Old mat file detected - skipping IGD/GD/Spread charts (no data)\n');
end

if isfield(results, 'wilcoxon_table_hv') && ~isempty(results.wilcoxon_table_hv)
    fprintf('正在生成 Wilcoxon 统计显著性表格 (HV + IGD + Spread)...\n');
    plotWilcoxonTableCombined(results, output_dir, scene);
end
end

function plotWilcoxonTableCombined(results, output_dir, scene)
    metric_names = {'HV', 'IGD', 'Spread'};
    metric_fields = {'wilcoxon_table_hv', 'wilcoxon_table_igd', 'wilcoxon_table_spread'};
    n_metrics = length(metric_names);

    tables = cell(n_metrics, 1);
    n_rows_all = zeros(n_metrics, 1);
    for m = 1:n_metrics
        if isfield(results, metric_fields{m}) && ~isempty(results.(metric_fields{m}))
            tables{m} = results.(metric_fields{m});
            n_rows_all(m) = size(tables{m}, 1);
        end
    end
    n_rows = max(n_rows_all);
    if n_rows == 0
        fprintf('  无 Wilcoxon 数据，跳过\n');
        return;
    end

    % --- 导出合并 CSV ---
    csv_file = fullfile(output_dir, 'wilcoxon_test_results.csv');
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Comparison');
    for m = 1:n_metrics
        fprintf(fid, ',%s_p-value,%s_significant', metric_names{m}, metric_names{m});
    end
    fprintf(fid, '\n');
    for r = 1:n_rows
        if ~isempty(tables{1}) && r <= size(tables{1}, 1)
            cmp_name = tables{1}{r, 1};
        elseif ~isempty(tables{2}) && r <= size(tables{2}, 1)
            cmp_name = tables{2}{r, 1};
        else
            cmp_name = '';
        end
        fprintf(fid, '%s', cmp_name);
        for m = 1:n_metrics
            if ~isempty(tables{m}) && r <= size(tables{m}, 1)
                p_val = tables{m}{r, 2};
                h_val = tables{m}{r, 3};
                if ~isnan(h_val) && h_val == 1
                    sig_str = 'Yes';
                elseif ~isnan(h_val)
                    sig_str = 'No';
                else
                    sig_str = 'N/A';
                end
                fprintf(fid, ',%.2e,%s', p_val, sig_str);
            else
                fprintf(fid, ',,');
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  Wilcoxon 检验结果已导出至: %s\n', csv_file);

    % --- 绘制合并表格图 ---
    fig_w = figure('Position', [100, 100, 950, 320 + n_rows * 32]);
    set(gcf, 'Color', 'w');

    col_names = {'Comparison', 'p-value', 'Sig.', 'p-value', 'Sig.', 'p-value', 'Sig.'};
    n_cols = length(col_names);
    table_data = cell(n_rows, n_cols);

    for r = 1:n_rows
        % 取对比名称，截短到 "cSA-GOA vs. XXX"
        if ~isempty(tables{1}) && r <= size(tables{1}, 1)
            raw_name = tables{1}{r, 1};
            paren_idx = find(raw_name == '(', 1, 'first');
            if ~isempty(paren_idx)
                raw_name = strtrim(raw_name(1:paren_idx-1));
            end
            table_data{r, 1} = raw_name;
        end
        col_idx = 2;
        for m = 1:n_metrics
            if ~isempty(tables{m}) && r <= size(tables{m}, 1)
                table_data{r, col_idx} = format_pvalue(tables{m}{r, 2});
                h_val = tables{m}{r, 3};
                if ~isnan(h_val) && h_val == 1
                    table_data{r, col_idx+1} = 'Yes';
                elseif ~isnan(h_val)
                    table_data{r, col_idx+1} = 'No';
                else
                    table_data{r, col_idx+1} = 'N/A';
                end
            else
                table_data{r, col_idx} = '-';
                table_data{r, col_idx+1} = '-';
            end
            col_idx = col_idx + 2;
        end
    end

    % 分组列颜色
    col_widths = {280, 110, 55, 110, 55, 110, 55};

    t = uitable('Data', table_data, 'ColumnName', col_names, ...
        'Position', [20, 40, 910, 240 + n_rows * 32], ...
        'FontName', 'Times New Roman', 'FontSize', 11, ...
        'ColumnWidth', col_widths, ...
        'RowName', [], ...
        'CellSelectionCallback', [], ...
        'ForegroundColor', [0.15, 0.15, 0.15]);

    % 着色：Sig列为Yes的标绿，No的标红
    % (uitable 不直接支持单元格着色，用 uitable 的 BackgroundColor 无法逐格设)
    % 改用 text-based 表格方案以获得更好的视觉效果

    close(fig_w);

    % --- 改用 axes + text 绘制，支持逐格着色 ---
    row_h = 28;
    header_h = 36;
    section_h = 28;
    fig_h = 80 + header_h + section_h + n_rows * row_h + 60;
    fig_w_px = 1100;
    fig_w2 = figure('Position', [100, 100, fig_w_px, fig_h], 'Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none');

    ax = axes('Position', [0, 0, 1, 1], 'XLim', [0, fig_w_px], 'YLim', [0, fig_h], ...
        'XTick', [], 'YTick', [], 'Color', 'w');
    hold(ax, 'on');

    % 列 x 坐标: Comparison | HV p | HV s | IGD p | IGD s | Spread p | Spread s
    col_x = [20, 420, 560, 640, 780, 860, 1000];
    col_align = {'left', 'center', 'center', 'center', 'center', 'center', 'center'};

    y_top = fig_h - 30;

    % 标题
    title_str = sprintf('Wilcoxon Rank Sum Test: cSA-GOA vs. Baselines  (%s)', scene);
    text(ax, fig_w_px/2, y_top, title_str, ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'center', 'Color', [0.1, 0.1, 0.1]);

    % 分组标题行 (HV / IGD / Spread)
    y_group = y_top - 30;
    group_centers = [15, (col_x(2)+col_x(3))/2, (col_x(4)+col_x(5))/2, (col_x(6)+col_x(7))/2];
    group_labels = {'Comparison', 'HV', 'IGD', 'Spread'};
    group_colors = {[0.3,0.3,0.3], [0.0, 0.35, 0.65], [0.65, 0.15, 0.15], [0.1, 0.5, 0.2]};

    % 画分组色条背景
    group_rects = {[col_x(1)-5, 400], [col_x(2)-10, col_x(3)+20], ...
                   [col_x(4)-10, col_x(5)+20], [col_x(6)-10, col_x(7)+60]};
    for g = 1:4
        if g > 1
            rectangle('Position', [group_rects{g}(1), y_group-12, group_rects{g}(2)-group_rects{g}(1), 24], ...
                'FaceColor', group_colors{g}*0.12 + 0.88*[1,1,1], 'EdgeColor', 'none', 'Curvature', 0.2);
        end
        text(ax, (group_rects{g}(1)+group_rects{g}(2))/2, y_group, group_labels{g}, ...
            'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'Color', group_colors{g});
    end

    % 子列标题
    y_header = y_group - section_h;
    sub_headers = {'Comparison', 'p-value', 'Sig.', 'p-value', 'Sig.', 'p-value', 'Sig.'};
    for c = 1:n_cols
        text(ax, col_x(c), y_header, sub_headers{c}, ...
            'FontSize', 10, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', col_align{c}, 'Color', [0.4, 0.4, 0.4]);
    end

    % 分隔线
    line(ax, [10, fig_w_px-10], [y_header-10, y_header-10], 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.5);

    % 数据行
    for r = 1:n_rows
        y_row = y_header - 18 - (r-1) * row_h;

        % 交替行背景
        if mod(r, 2) == 0
            rectangle('Position', [5, y_row-10, fig_w_px-10, row_h], ...
                'FaceColor', [0.96, 0.96, 0.98], 'EdgeColor', 'none');
        end

        % 对比名称
        text(ax, col_x(1), y_row, table_data{r, 1}, ...
            'FontSize', 10.5, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'left', 'Color', [0.15, 0.15, 0.15]);

        col_idx = 2;
        for m = 1:n_metrics
            % p-value
            text(ax, col_x(col_idx), y_row, table_data{r, col_idx}, ...
                'FontSize', 10, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'Color', [0.2, 0.2, 0.2]);

            % Sig.
            sig_text = table_data{r, col_idx+1};
            if strcmp(sig_text, 'Yes')
                sig_color = [0.0, 0.55, 0.0];  % 绿色
                sig_display = 'Yes';
            elseif strcmp(sig_text, 'No')
                sig_color = [0.7, 0.15, 0.1];  % 红色
                sig_display = 'No';
            else
                sig_color = [0.5, 0.5, 0.5];
                sig_display = sig_text;
            end
            text(ax, col_x(col_idx+1), y_row, sig_display, ...
                'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'Color', sig_color);

            col_idx = col_idx + 2;
        end

        % 行分隔线
        line(ax, [10, fig_w_px-10], [y_row-12, y_row-12], 'Color', [0.88, 0.88, 0.88], 'LineWidth', 0.5);
    end

    % 底部注释
    y_note = y_header - 25 - n_rows * row_h;
    text(ax, 15, y_note, 'Sig.: significant at \alpha = 0.05  |  Yes = cSA-GOA significantly different (and better for + metrics)', ...
        'FontSize', 9, 'FontName', 'Times New Roman', 'Color', [0.5, 0.5, 0.5], ...
        'HorizontalAlignment', 'left');

    hold(ax, 'off');
    axis(ax, 'off');

    saveas(fig_w2, fullfile(output_dir, 'wilcoxon_table.fig'));
    exportgraphics(fig_w2, fullfile(output_dir, 'wilcoxon_table.png'), 'Resolution', 300);
    close(fig_w2);
    fprintf('  Wilcoxon 统计表格已保存至: %s\n', output_dir);
end

function s = format_pvalue(p)
    if isnan(p)
        s = 'N/A';
    elseif p == 0
        s = '< 1.00e-300';
    elseif p < 0.001
        s = sprintf('%.2e', p);
    else
        s = sprintf('%.6f', p);
    end
end
