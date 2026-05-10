function plot_ablation_all(data_file, output_dir)
% plot_ablation_all - 娑堣瀺瀹為獙瀹屾暣鍙鍖栵紙鏀舵暃鏇茬嚎/HV鏌辩姸鍥?绠卞瀷鍥?缁熻琛級
% 支持自动识别场景（Small/Medium/Large）并输出到对应子文件夹
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'ablation_results_para_*.mat');
    files = dir(data_file);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        data_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No ablation results file found');
    end
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures', 'ablation');
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

variant_order = {'proposed', 'no_pareto', 'no_subpop', 'no_goa', 'no_migration'};
variant_labels_map = struct( ...
    'proposed', 'Proposed cSA-GOA', ...
    'no_pareto', 'w/o Pareto Leader', ...
    'no_subpop', 'w/o Multi-Subpop', ...
    'no_goa', 'w/o GOA Repulsion', ...
    'no_migration', 'w/o Elite Migration');
variant_colors_map = struct( ...
    'proposed', [0.8500, 0.3250, 0.0980], ...
    'no_pareto', [0.0000, 0.4470, 0.7410], ...
    'no_subpop', [0.9290, 0.6940, 0.1250], ...
    'no_goa', [0.4940, 0.1840, 0.5560], ...
    'no_migration', [0.4660, 0.6740, 0.1880]);
variant_markers_map = struct( ...
    'proposed', 'o', ...
    'no_pareto', 's', ...
    'no_subpop', '^', ...
    'no_goa', 'd', ...
    'no_migration', 'p');
variant_size_map = struct( ...
    'proposed', 30, ...
    'no_pareto', 20, ...
    'no_subpop', 20, ...
    'no_goa', 20, ...
    'no_migration', 20);

variants = {};
labels = {};
colors = zeros(0, 3);
marker_styles = {};
marker_sizes = [];
for idx = 1:length(variant_order)
    variant_name = variant_order{idx};
    if isfield(results, variant_name)
        variants{end + 1} = variant_name; %#ok<AGROW>
        labels{end + 1} = variant_labels_map.(variant_name); %#ok<AGROW>
        colors(end + 1, :) = variant_colors_map.(variant_name); %#ok<AGROW>
        marker_styles{end + 1} = variant_markers_map.(variant_name); %#ok<AGROW>
        marker_sizes(end + 1) = variant_size_map.(variant_name); %#ok<AGROW>
    end
end
n_vars = length(variants);

%% ========== 鍥?: 鏀舵暃鏇茬嚎 + HV鏌辩姸鍥?==========
fig1 = figure('Position', [100, 100, 900, 500]);
set(gcf, 'Color', 'w');

subplot(1, 2, 1);
hold on;
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'convergence_curves')
        curves = results.(variant).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            max_len = max(cellfun(@length, curves));
            avg_curve = zeros(1, max_len);
            n_runs = 0;
            for run = 1:min(length(curves), 30)
                if ~isempty(curves{run})
                    c = curves{run};
                    avg_curve(1:length(c)) = avg_curve(1:length(c)) + c;
                    n_runs = n_runs + 1;
                end
            end
            if n_runs > 0
                avg_curve = avg_curve / n_runs;
                lw = 2.0;
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors(v,:), 'LineWidth', lw, 'DisplayName', labels{v});
            end
        end
    end
end
xlabel('Generations', 'FontWeight', 'bold');
ylabel('Best Fitness', 'FontWeight', 'bold');
title(['Convergence Curves (Average of 30 runs) - ', scene], 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; box on;

subplot(1, 2, 2);
hv_means = zeros(1, n_vars);
hv_stds = zeros(1, n_vars);
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'hv_values')
        vals = results.(variant).hv_values;
        hv_means(v) = nanmean(vals);
        hv_stds(v) = nanstd(vals);
    end
end
% 鉁?淇鐐?2锛氬皢鍐欐鐨?5 鍏ㄩ儴鏇挎崲涓哄姩鎬佺殑 n_vars
bar(1:n_vars, hv_means, 0.6, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'k', 'LineWidth', 1.2);
hold on;
errorbar(1:n_vars, hv_means, hv_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'Capsize', 6);
xlabel('Variant', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title(['HV Comparison (Mean 卤 Std) - ', scene], 'FontWeight', 'bold');
xticks(1:n_vars);
xticklabels(labels);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; box on;

saveas(fig1, fullfile(output_dir, 'ablation_convergence_hv.fig'));
saveas(fig1, fullfile(output_dir, 'ablation_convergence_hv.png'));
close(fig1);

%% ========== 鍥?: 鍙孻杞存煴鐘跺浘 (鎬ц兘 vs 璁＄畻鎴愭湰) ==========
fig2 = figure('Position', [100, 100, 800, 500]);
set(gcf, 'Color', 'w');

hv_means = zeros(1, n_vars);
iter_means = zeros(1, n_vars);
for i = 1:n_vars
    v = variants{i};
    if isfield(results, v)
        hv_means(i) = results.(v).mean_hv;
        iter_means(i) = mean(results.(v).iter_counts);
    end
end

yyaxis left
b1 = bar(1:n_vars, hv_means, 0.4, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold', 'Color', [0.0000, 0.4470, 0.7410]);
ylim([min(hv_means)*0.95, max(hv_means)*1.05]);
set(gca, 'ycolor', [0.0000, 0.4470, 0.7410]);

yyaxis right
hold on;
b2 = bar((1:n_vars)+0.4, iter_means, 0.4, 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Actual Iterations (Cost)', 'FontWeight', 'bold', 'Color', [0.8500, 0.3250, 0.0980]);
ylim([0, 350]);
set(gca, 'ycolor', [0.8500, 0.3250, 0.0980]);

xticks((1:n_vars) + 0.2);
xticklabels(labels);
legend([b1, b2], {'Performance (HV)', 'Computational Cost (Iter)'}, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on; box on;

saveas(fig2, fullfile(output_dir, 'ablation_performance_cost.fig'));
saveas(fig2, fullfile(output_dir, 'ablation_performance_cost.png'));
close(fig2);

%% ========== 鍥?: 绠卞瀷鍥?(HV鍒嗗竷) ==========
fig3 = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

all_hv = [];
group_idx = [];
for i = 1:n_vars
    v = variants{i};
    if isfield(results, v) && isfield(results.(v), 'hv_values')
        hv_data = results.(v).hv_values;
        all_hv = [all_hv; hv_data(:)];
        group_idx = [group_idx; i * ones(length(hv_data(:)), 1)];
    end
end

hold on;
for i = 1:n_vars
    idx = (group_idx == i);
    if any(idx)
        b = boxchart(group_idx(idx), all_hv(idx));
        b.BoxFaceColor = colors(i, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:n_vars);
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title(['HV Distribution (Robustness Analysis - 30 Runs) - ', scene], 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3; box on;

saveas(fig3, fullfile(output_dir, 'ablation_boxplot.fig'));
saveas(fig3, fullfile(output_dir, 'ablation_boxplot.png'));
close(fig3);

%% ========== 鍥?: 缁熻妫€楠岃〃 ==========
fig4 = figure('Position', [100, 100, 700, 300]);
set(gcf, 'Color', 'w');

variant_labels = labels;

if ~isfield(results, 'proposed')
    fprintf('Error: proposed variant not found in data\n');
    close(fig4);
    return;
end

proposed_hv = results.proposed.hv_values;

p_values = zeros(n_vars, 1);
h_stats = zeros(n_vars, 1);
mean_improvement = zeros(n_vars, 1);

for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && strcmp(variant, 'proposed')
        p_values(v) = 1.0;
        h_stats(v) = 0;
        mean_improvement(v) = 0;
    elseif isfield(results, variant)
        variant_hv = results.(variant).hv_values;
        valid_idx = ~isnan(proposed_hv) & ~isnan(variant_hv);
        if sum(valid_idx) >= 3
            [p_values(v), h_stats(v)] = local_signrank(proposed_hv(valid_idx), variant_hv(valid_idx));
        else
            p_values(v) = NaN;
            h_stats(v) = NaN;
        end
        mean_improvement(v) = (nanmean(proposed_hv) - nanmean(variant_hv)) / nanmean(variant_hv) * 100;
    else
        p_values(v) = NaN; h_stats(v) = NaN; mean_improvement(v) = NaN;
    end
end

sig_labels = cell(n_vars, 1);
for v = 1:n_vars
    if p_values(v) < 0.001
        sig_labels{v} = '***';
    elseif p_values(v) < 0.01
        sig_labels{v} = '**';
    elseif p_values(v) < 0.05
        sig_labels{v} = '*';
    else
        sig_labels{v} = 'ns';
    end
end

col_labels = {'Variant', 'HV (Mean +/- Std)', 'p-value', 'Sig.', 'Improve (%)'};
table_data = cell(n_vars, 5);
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant)
        hv_mean = nanmean(results.(variant).hv_values);
        hv_std = nanstd(results.(variant).hv_values);
        table_data{v, 1} = variant_labels{v};
        table_data{v, 2} = sprintf('%.4f +/- %.4f', hv_mean, hv_std);
        table_data{v, 3} = sprintf('%.4e', p_values(v));
        table_data{v, 4} = sig_labels{v};
        table_data{v, 5} = sprintf('%.2f%%', mean_improvement(v));
    end
end

axis off;
title(['Statistical Tests (Wilcoxon Signed-Rank Test, 30 Runs) - ', scene], ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 13);

n_cols = 5;
n_rows = n_vars + 1;
col_widths = [0.28, 0.22, 0.18, 0.10, 0.16];
x_start = 0.03;
y_top = 0.88;
row_height = 0.85 / n_rows;

for c = 1:n_cols
    x_pos = x_start + sum(col_widths(1:c-1));
    text(x_pos, y_top, col_labels{c}, 'Units', 'normalized', ...
        'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

line_x = [x_start, x_start + sum(col_widths)];
line_y = y_top - row_height * 0.5;
annotation('line', line_x, [line_y, line_y], 'LineWidth', 1.5, 'Color', [0.3, 0.3, 0.3]);

for r = 1:n_vars
    y_pos = y_top - r * row_height;
    for c = 1:n_cols
        x_pos = x_start + sum(col_widths(1:c-1));
        cell_text = table_data{r, c};
        if c == 4
            if strcmp(cell_text, '***')
                txt_color = [0.8, 0, 0];
            elseif strcmp(cell_text, '**')
                txt_color = [0.9, 0.4, 0];
            elseif strcmp(cell_text, '*')
                txt_color = [0, 0.5, 0];
            else
                txt_color = [0.4, 0.4, 0.4];
            end
        else
            txt_color = [0, 0, 0];
        end
        text(x_pos, y_pos, cell_text, 'Units', 'normalized', ...
            'FontName', 'Times New Roman', 'FontSize', 10, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'Color', txt_color);
    end
    if r < n_vars
        sep_y = y_pos - row_height * 0.5;
        annotation('line', line_x, [sep_y, sep_y], 'LineWidth', 0.5, ...
            'LineStyle', ':', 'Color', [0.7, 0.7, 0.7]);
    end
end

exportgraphics(fig4, fullfile(output_dir, 'ablation_statistical_tests.png'), 'Resolution', 300);
saveas(fig4, fullfile(output_dir, 'ablation_statistical_tests.fig'));
close(fig4);

fig5 = figure('Position', [100, 100, 800, 600]);
set(gcf, 'Color', 'w');
view(45, 30);
hold on;

sizes = marker_sizes + 20;

for v = 1:n_vars
    variant = variants{v};
    if ~isfield(results, variant) || ~isfield(results.(variant), 'pareto_fronts')
        continue;
    end
    pfs = results.(variant).pareto_fronts;
    if isempty(pfs)
        continue;
    end
    all_u = []; all_l = []; all_e = [];
    for run = 1:length(pfs)
        pf = pfs{run};
        if ~isempty(pf) && size(pf, 2) == 3
            all_u = [all_u; pf(:,1)];
            all_l = [all_l; pf(:,2)];
            all_e = [all_e; pf(:,3)];
        end
    end

    if ~isempty(all_u)
        scatter3(all_l, all_e, all_u, sizes(v), marker_styles{v}, ...
            'filled', 'MarkerFaceColor', colors(v,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
            'DisplayName', labels{v});
    end
end

xlabel('Total Latency (s) \downarrow', 'FontWeight', 'bold');
ylabel('Total Energy (J) \downarrow', 'FontWeight', 'bold');
zlabel('System Utility \uparrow', 'FontWeight', 'bold');
title(['3D Pareto Front - Ablation Study (', scene, ')'], 'FontWeight', 'bold');
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; box on;
view(45, 30);

saveas(fig5, fullfile(output_dir, 'ablation_pareto_3d.fig'));
saveas(fig5, fullfile(output_dir, 'ablation_pareto_3d.png'));
close(fig5);

fig6 = figure('Name', 'Pareto Projections', 'Position', [100, 100, 1500, 450]);
set(gcf, 'Color', 'w');

proj_configs = {
    2, 1, 'Total Latency (s) \downarrow', 'System Utility \uparrow', 1;
    3, 1, 'Total Energy (J) \downarrow', 'System Utility \uparrow', 2;
    2, 3, 'Total Latency (s) \downarrow', 'Total Energy (J) \downarrow', 3
};

for p = 1:3
    subplot(1, 3, proj_configs{p, 5});
    hold on; grid on;

    for v = n_vars:-1:1
        variant = variants{v};
        if ~isfield(results, variant) || ~isfield(results.(variant), 'pareto_fronts'), continue; end

        pfs = results.(variant).pareto_fronts;
        all_data = [];
        for run = 1:length(pfs)
            if ~isempty(pfs{run}), all_data = [all_data; pfs{run}]; end
        end

        if ~isempty(all_data)
            if strcmp(variant, 'proposed')
                h = scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    marker_sizes(v), marker_styles{v}, ...
                    'MarkerEdgeColor', colors(v,:), ...
                    'LineWidth', 1.1, ...
                    'DisplayName', labels{v});
                if isprop(h, 'MarkerEdgeAlpha')
                    h.MarkerEdgeAlpha = 0.58;
                end
            else
                h = scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    marker_sizes(v), colors(v,:), marker_styles{v}, 'filled', ...
                    'MarkerFaceAlpha', 0.42, ...
                    'MarkerEdgeColor', colors(v,:), ...
                    'LineWidth', 0.35, ...
                    'DisplayName', labels{v});
                if isprop(h, 'MarkerEdgeAlpha')
                    h.MarkerEdgeAlpha = 0.28;
                end
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

saveas(fig6, fullfile(output_dir, 'ablation_projections_2d.png'));
saveas(fig6, fullfile(output_dir, 'ablation_projections_2d.fig'));
close(fig6);

fprintf('Ablation visualization saved to %s\n', output_dir);
end

function [p, h] = local_signrank(x, y)
    try
        [p, h] = signrank(x, y);
    catch
        [p, h] = wilcoxon_signed_rank_selfcontained(x, y);
    end
end

function [p, h] = wilcoxon_signed_rank_selfcontained(x, y)
    d = x(:) - y(:);
    d = d(d ~= 0);
    n = length(d);
    if n < 3
        p = NaN; h = NaN;
        return;
    end
    [~, idx] = sort(abs(d));
    ranks = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && abs(d(idx(j+1))) == abs(d(idx(i)))
            j = j + 1;
        end
        ranks(idx(i:j)) = (i + j) / 2;
        i = j + 1;
    end
    W = sum(sign(d) .* ranks);
    sigma_w = sqrt(n * (n + 1) * (2 * n + 1) / 6);
    if sigma_w > 0
        z = (abs(W) - 0.5) / sigma_w;
        p = 2 * (1 - normcdf(z));
    else
        p = 1;
    end
    h = double(p < 0.05);
end

