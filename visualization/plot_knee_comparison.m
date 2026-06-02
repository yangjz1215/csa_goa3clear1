function plot_knee_comparison(data_file, output_dir)
% plot_knee_comparison - 膝点解提取与多维度对比可视化
% 包含：(1) 膝点三目标表格  (2) 归一化综合得分  (3) 极端解对比  (4) 雷达图
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

algorithms = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA2'};
labels = {'cSA-GOA (Proposed)', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA-II'};
n_algs = length(algorithms);

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

knee_utils = nan(n_algs, 30);
knee_lats  = nan(n_algs, 30);
knee_nrgs  = nan(n_algs, 30);

best_util_utils = nan(n_algs, 30);
best_util_lats  = nan(n_algs, 30);
best_util_nrgs  = nan(n_algs, 30);

best_lat_utils = nan(n_algs, 30);
best_lat_lats  = nan(n_algs, 30);
best_lat_nrgs  = nan(n_algs, 30);

best_nrg_utils = nan(n_algs, 30);
best_nrg_lats  = nan(n_algs, 30);
best_nrg_nrgs  = nan(n_algs, 30);

for a = 1:n_algs
    alg = algorithms{a};
    if ~isfield(results, alg) || ~isfield(results.(alg), 'pareto_fronts')
        continue;
    end
    pfs = results.(alg).pareto_fronts;
    n_runs = length(pfs);
    for run = 1:n_runs
        pf = pfs{run};
        if isempty(pf) || size(pf, 1) < 1
            continue;
        end

        if size(pf, 1) >= 2
            knee_idx = findKneePoint(pf);
        else
            knee_idx = 1;
        end
        knee_utils(a, run) = pf(knee_idx, 1);
        knee_lats(a, run)  = pf(knee_idx, 2);
        knee_nrgs(a, run)  = pf(knee_idx, 3);

        [~, idx_u] = max(pf(:, 1));
        best_util_utils(a, run) = pf(idx_u, 1);
        best_util_lats(a, run)  = pf(idx_u, 2);
        best_util_nrgs(a, run)  = pf(idx_u, 3);

        [~, idx_l] = min(pf(:, 2));
        best_lat_utils(a, run) = pf(idx_l, 1);
        best_lat_lats(a, run)  = pf(idx_l, 2);
        best_lat_nrgs(a, run)  = pf(idx_l, 3);

        [~, idx_e] = min(pf(:, 3));
        best_nrg_utils(a, run) = pf(idx_e, 1);
        best_nrg_lats(a, run)  = pf(idx_e, 2);
        best_nrg_nrgs(a, run)  = pf(idx_e, 3);
    end
end

mean_ku = mean(knee_utils, 2, 'omitnan');
std_ku  = std(knee_utils, 0, 2, 'omitnan');
mean_kl = mean(knee_lats, 2, 'omitnan');
std_kl  = std(knee_lats, 0, 2, 'omitnan');
mean_ke = mean(knee_nrgs, 2, 'omitnan');
std_ke  = std(knee_nrgs, 0, 2, 'omitnan');

all_utils = [knee_utils; best_util_utils; best_lat_utils; best_nrg_utils];
all_lats  = [knee_lats;  best_util_lats;  best_lat_lats;  best_nrg_lats];
all_nrgs  = [knee_nrgs;  best_util_nrgs;  best_lat_nrgs;  best_nrg_nrgs];

global_max_u = max(all_utils(:), [], 'omitnan');
global_min_u = min(all_utils(:), [], 'omitnan');
global_max_l = max(all_lats(:), [], 'omitnan');
global_min_l = min(all_lats(:), [], 'omitnan');
global_max_e = max(all_nrgs(:), [], 'omitnan');
global_min_e = min(all_nrgs(:), [], 'omitnan');

range_u = global_max_u - global_min_u; if range_u < 1e-10, range_u = 1; end
range_l = global_max_l - global_min_l; if range_l < 1e-10, range_l = 1; end
range_e = global_max_e - global_min_e; if range_e < 1e-10, range_e = 1; end

norm_ku = (mean_ku - global_min_u) / range_u;
norm_kl = (global_max_l - mean_kl) / range_l;
norm_ke = (global_max_e - mean_ke) / range_e;

composite_score = sqrt((1 - norm_ku).^2 + (1 - norm_kl).^2 + (1 - norm_ke).^2);

energy_efficiency = mean_ku ./ max(mean_ke, 1e-10);

fprintf('\n========== Knee Point Comparison (%s) ==========\n', scene);
fprintf('%-22s | %-24s | %-24s | %-24s | %-12s | %-14s\n', 'Algorithm', 'Utility (↑)', 'Latency/s (↓)', 'Energy/J (↓)', 'D_ideal (↓)', 'U/E (↑)');
fprintf('%s\n', repmat('-', 1, 140));
for a = 1:n_algs
    fprintf('%-22s | %8.2f ± %-8.2f | %8.4f ± %-8.4f | %8.2f ± %-8.2f | %8.4f | %8.6f\n', ...
        labels{a}, mean_ku(a), std_ku(a), mean_kl(a), std_kl(a), mean_ke(a), std_ke(a), composite_score(a), energy_efficiency(a));
end
fprintf('%s\n', repmat('-', 1, 140));

[~, best_idx] = min(composite_score);
fprintf('  >> Best knee point (closest to ideal): %s (D_ideal = %.4f)\n', labels{best_idx}, composite_score(best_idx));
[~, best_eff_idx] = max(energy_efficiency);
fprintf('  >> Best energy efficiency: %s (U/E = %.6f)\n', labels{best_eff_idx}, energy_efficiency(best_eff_idx));

csv_file = fullfile(output_dir, 'knee_point_comparison.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'Algorithm,Utility_Mean,Utility_Std,Latency_Mean,Latency_Std,Energy_Mean,Energy_Std,Norm_Utility,Norm_Latency,Norm_Energy,Composite_D_Ideal,Energy_Efficiency\n');
for a = 1:n_algs
    fprintf(fid, '%s,%.4f,%.4f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.8f\n', ...
        labels{a}, mean_ku(a), std_ku(a), mean_kl(a), std_kl(a), mean_ke(a), std_ke(a), ...
        norm_ku(a), norm_kl(a), norm_ke(a), composite_score(a), energy_efficiency(a));
end
fclose(fid);
fprintf('Knee point CSV saved to: %s\n', csv_file);

drawKneeTable(labels, mean_ku, std_ku, mean_kl, std_kl, mean_ke, std_ke, composite_score, energy_efficiency, colors, output_dir, scene);

drawKneeBarChart(labels, mean_ku, std_ku, mean_kl, std_kl, mean_ke, std_ke, colors, output_dir, scene);

drawRadarChart(labels, norm_ku, norm_kl, norm_ke, composite_score, colors, output_dir, scene);

drawEfficiencyRadar(labels, norm_ku, norm_kl, norm_ke, energy_efficiency, composite_score, colors, output_dir, scene);

drawExtremeComparison(labels, ...
    mean(best_util_utils, 2, 'omitnan'), std(best_util_utils, 0, 2, 'omitnan'), ...
    mean(best_util_lats, 2, 'omitnan'),  std(best_util_lats, 0, 2, 'omitnan'), ...
    mean(best_util_nrgs, 2, 'omitnan'),  std(best_util_nrgs, 0, 2, 'omitnan'), ...
    mean(best_lat_utils, 2, 'omitnan'),  std(best_lat_utils, 0, 2, 'omitnan'), ...
    mean(best_lat_lats, 2, 'omitnan'),   std(best_lat_lats, 0, 2, 'omitnan'), ...
    mean(best_lat_nrgs, 2, 'omitnan'),   std(best_lat_nrgs, 0, 2, 'omitnan'), ...
    mean(best_nrg_utils, 2, 'omitnan'),  std(best_nrg_utils, 0, 2, 'omitnan'), ...
    mean(best_nrg_lats, 2, 'omitnan'),   std(best_nrg_lats, 0, 2, 'omitnan'), ...
    mean(best_nrg_nrgs, 2, 'omitnan'),   std(best_nrg_nrgs, 0, 2, 'omitnan'), ...
    colors, output_dir, scene);

fprintf('Knee point comparison charts saved to: %s\n', output_dir);
end

function knee_idx = findKneePoint(pf)
    n = size(pf, 1);
    if n == 1
        knee_idx = 1;
        return;
    end

    max_u = max(pf(:, 1)); min_u = min(pf(:, 1));
    max_l = max(pf(:, 2)); min_l = min(pf(:, 2));
    max_e = max(pf(:, 3)); min_e = min(pf(:, 3));

    range_u = max_u - min_u; if range_u < 1e-10, range_u = 1; end
    range_l = max_l - min_l; if range_l < 1e-10, range_l = 1; end
    range_e = max_e - min_e; if range_e < 1e-10, range_e = 1; end

    pf_norm = zeros(n, 3);
    pf_norm(:, 1) = (pf(:, 1) - min_u) / range_u;
    pf_norm(:, 2) = (max_l - pf(:, 2)) / range_l;
    pf_norm(:, 3) = (max_e - pf(:, 3)) / range_e;

    ideal = [1, 1, 1];
    w = [1/3, 1/3, 1/3];
    rho = 0.001;

    cheb_dist = zeros(n, 1);
    for i = 1:n
        gap = w .* abs(ideal - pf_norm(i, :));
        cheb_dist(i) = max(gap) + rho * sum(gap);
    end

    [~, knee_idx] = min(cheb_dist);
end

function drawKneeTable(labels, mean_ku, std_ku, mean_kl, std_kl, mean_ke, std_ke, composite_score, energy_efficiency, colors, output_dir, scene)
    n_algs = length(labels);
    row_h = 32;
    header_h = 40;
    title_h = 50;
    fig_w_px = 1400;
    fig_h = title_h + header_h + n_algs * row_h + 40;

    fig = figure('Position', [100, 100, fig_w_px, fig_h], 'Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Position', [0, 0, 1, 1], 'XLim', [0, fig_w_px], 'YLim', [0, fig_h], ...
        'XTick', [], 'YTick', [], 'Color', 'w');
    hold(ax, 'on');

    col_x = [30, 260, 470, 680, 890, 1090];
    col_labels = {'Algorithm', 'Utility (\uparrow)', 'Latency / s (\downarrow)', 'Energy / J (\downarrow)', 'D_{ideal} (\downarrow)', 'U/E (\uparrow)'};
    col_widths_px = [230, 210, 210, 210, 200, 200];

    y_cursor = fig_h - 15;
    text(ax, fig_w_px / 2, y_cursor, sprintf('Knee Point Solution Comparison (%s Scene)', scene), ...
        'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'center', 'Color', [0.1, 0.1, 0.1]);

    y_header = y_cursor - header_h + 8;

    header_rect_y = y_header - 14;
    rectangle('Position', [col_x(1) - 10, header_rect_y, sum(col_widths_px) + 20, 28], ...
        'FaceColor', [0.15, 0.30, 0.55], 'EdgeColor', [0.1, 0.2, 0.4], 'Curvature', 0.1);

    for c = 1:6
        x_center = col_x(c) + col_widths_px(c) / 2;
        text(ax, x_center, y_header, col_labels{c}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [1, 1, 1]);
    end

    [~, best_comp_idx] = min(composite_score);
    [~, best_eff_idx] = max(energy_efficiency);

    for a = 1:n_algs
        y_row = y_header - 14 - (a - 0.5) * row_h;

        if a == best_comp_idx
            row_bg_color = [0.92, 0.98, 0.92];
        elseif mod(a, 2) == 0
            row_bg_color = [0.95, 0.95, 0.97];
        else
            row_bg_color = [1, 1, 1];
        end
        rectangle('Position', [col_x(1) - 10, y_row - row_h/2, sum(col_widths_px) + 20, row_h], ...
            'FaceColor', row_bg_color, 'EdgeColor', [0.85, 0.85, 0.85], 'Curvature', 0);

        text(ax, col_x(1) + 8, y_row, labels{a}, ...
            'FontSize', 10.5, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'Color', colors(a, :));

        util_str = sprintf('%.2f \\pm %.2f', mean_ku(a), std_ku(a));
        lat_str  = sprintf('%.4f \\pm %.4f', mean_kl(a), std_kl(a));
        nrg_str  = sprintf('%.2f \\pm %.2f', mean_ke(a), std_ke(a));
        comp_str = sprintf('%.4f', composite_score(a));
        eff_str  = sprintf('%.6f', energy_efficiency(a));

        text(ax, col_x(2) + col_widths_px(2)/2, y_row, util_str, ...
            'FontSize', 10.5, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', [0.15, 0.15, 0.15]);

        text(ax, col_x(3) + col_widths_px(3)/2, y_row, lat_str, ...
            'FontSize', 10.5, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', [0.15, 0.15, 0.15]);

        text(ax, col_x(4) + col_widths_px(4)/2, y_row, nrg_str, ...
            'FontSize', 10.5, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', [0.15, 0.15, 0.15]);

        comp_color = [0.15, 0.15, 0.15];
        if a == best_comp_idx
            comp_color = [0.0, 0.5, 0.0];
        end
        text(ax, col_x(5) + col_widths_px(5)/2, y_row, comp_str, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', comp_color);

        eff_color = [0.15, 0.15, 0.15];
        if a == best_eff_idx
            eff_color = [0.0, 0.5, 0.0];
        end
        text(ax, col_x(6) + col_widths_px(6)/2, y_row, eff_str, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', eff_color);
    end

    saveas(fig, fullfile(output_dir, 'knee_point_table.fig'));
    saveas(fig, fullfile(output_dir, 'knee_point_table.png'));
    close(fig);
end

function drawKneeBarChart(labels, mean_ku, std_ku, mean_kl, std_kl, mean_ke, std_ke, colors, output_dir, scene)
    n_algs = length(labels);

    fig = figure('Position', [100, 100, 1400, 420], 'Color', 'w');

    obj_titles = {'System Utility (\uparrow)', 'Total Latency / s (\downarrow)', 'Total Energy / J (\downarrow)'};
    means = {mean_ku, mean_kl, mean_ke};
    stds  = {std_ku, std_kl, std_ke};

    for obj = 1:3
        subplot(1, 3, obj);
        hold on;

        m = means{obj};
        s = stds{obj};

        for a = 1:n_algs
            bar(a, m(a), 0.6, 'FaceColor', colors(a, :), 'EdgeColor', 'k', 'LineWidth', 1.0, 'FaceAlpha', 0.85);
        end
        errorbar(1:n_algs, m, s, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);

        set(gca, 'XTick', 1:n_algs, 'XTickLabel', labels, 'FontName', 'Times New Roman', 'FontSize', 9);
        xtickangle(30);
        ylabel(obj_titles{obj}, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 11);
        title(['Knee Point: ', obj_titles{obj}], 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
        grid on;
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
        box on;
    end

    sgtitle(sprintf('Knee Point Solution Comparison (%s Scene)', scene), ...
        'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 14);

    saveas(fig, fullfile(output_dir, 'knee_point_bar_chart.fig'));
    saveas(fig, fullfile(output_dir, 'knee_point_bar_chart.png'));
    close(fig);
end

function drawRadarChart(labels, norm_ku, norm_kl, norm_ke, composite_score, colors, output_dir, scene)
    n_algs = length(labels);

    fig = figure('Position', [100, 100, 700, 600], 'Color', 'w');
    hold on;

    n_axes = 3;
    angles = linspace(0, 2*pi, n_axes + 1);
    angles = angles(1:end-1);

    axis_labels = {'Utility', 'Latency', 'Energy'};
    axis_labels_full = {'Utility (\uparrow)', 'Latency (\downarrow)', 'Energy (\downarrow)'};

    for ring = [0.2, 0.4, 0.6, 0.8, 1.0]
        ring_x = ring * cos(angles);
        ring_y = ring * sin(angles);
        ring_x = [ring_x, ring_x(1)];
        ring_y = [ring_y, ring_y(1)];
        plot(ring_x, ring_y, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);
    end

    for k = 1:n_axes
        plot([0, cos(angles(k))], [0, sin(angles(k))], '-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.8);
    end

    for a = 1:n_algs
        vals = [norm_ku(a), norm_kl(a), norm_ke(a)];
        x_pts = vals .* cos(angles);
        y_pts = vals .* sin(angles);
        x_pts = [x_pts, x_pts(1)];
        y_pts = [y_pts, y_pts(1)];

        if a == 1
            lw = 2.5;
            fa = 0.15;
        else
            lw = 1.2;
            fa = 0.05;
        end

        fill(x_pts, y_pts, colors(a, :), 'FaceAlpha', fa, 'EdgeColor', colors(a, :), 'LineWidth', lw, ...
            'DisplayName', sprintf('%s (D=%.3f)', labels{a}, composite_score(a)));
    end

    for k = 1:n_axes
        label_r = 1.18;
        text(label_r * cos(angles(k)), label_r * sin(angles(k)), axis_labels_full{k}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    axis equal;
    axis off;
    title(sprintf('Normalized Knee Point Radar (%s Scene)', scene), ...
        'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 14);
    legend('Location', 'southoutside', 'FontName', 'Times New Roman', 'FontSize', 9, 'NumColumns', 2);

    saveas(fig, fullfile(output_dir, 'knee_point_radar.fig'));
    saveas(fig, fullfile(output_dir, 'knee_point_radar.png'));
    close(fig);
end

function drawExtremeComparison(labels, ...
    bu_u, bu_u_s, bu_l, bu_l_s, bu_e, bu_e_s, ...
    bl_u, bl_u_s, bl_l, bl_l_s, bl_e, bl_e_s, ...
    be_u, be_u_s, be_l, be_l_s, be_e, be_e_s, ...
    colors, output_dir, scene)

    n_algs = length(labels);
    row_h = 28;
    header_h = 36;
    section_h = 30;
    title_h = 50;
    fig_w_px = 1200;

    n_rows = n_algs * 3;
    fig_h = title_h + section_h * 3 + header_h * 3 + n_rows * row_h + 60;

    fig = figure('Position', [50, 50, fig_w_px, min(fig_h, 900)], 'Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Position', [0, 0, 1, 1], 'XLim', [0, fig_w_px], 'YLim', [0, fig_h], ...
        'XTick', [], 'YTick', [], 'Color', 'w');
    hold(ax, 'on');

    col_x = [30, 280, 530, 780];
    col_widths_px = [250, 250, 250, 250];
    col_labels = {'Algorithm', 'Utility (\uparrow)', 'Latency / s (\downarrow)', 'Energy / J (\downarrow)'};

    y_cursor = fig_h - 20;
    text(ax, fig_w_px / 2, y_cursor, sprintf('Extreme Solution Comparison (%s Scene)', scene), ...
        'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'center', 'Color', [0.1, 0.1, 0.1]);

    sections = {
        'Best Utility Solution (max Utility from Pareto front)', bu_u, bu_u_s, bu_l, bu_l_s, bu_e, bu_e_s;
        'Best Latency Solution (min Latency from Pareto front)', bl_u, bl_u_s, bl_l, bl_l_s, bl_e, bl_e_s;
        'Best Energy Solution (min Energy from Pareto front)', be_u, be_u_s, be_l, be_l_s, be_e, be_e_s
    };

    for sec = 1:3
        y_cursor = y_cursor - section_h;
        sec_color = [0.15, 0.30, 0.55] * (0.7 + 0.1 * sec);
        rectangle('Position', [col_x(1) - 10, y_cursor - 12, sum(col_widths_px) + 20, 24], ...
            'FaceColor', sec_color, 'EdgeColor', 'none', 'Curvature', 0.1);
        text(ax, fig_w_px / 2, y_cursor, sections{sec, 1}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [1, 1, 1]);

        y_cursor = y_cursor - header_h + 14;
        rectangle('Position', [col_x(1) - 10, y_cursor - 10, sum(col_widths_px) + 20, 20], ...
            'FaceColor', [0.85, 0.88, 0.92], 'EdgeColor', 'none');
        for c = 1:4
            x_center = col_x(c) + col_widths_px(c) / 2;
            text(ax, x_center, y_cursor, col_labels{c}, ...
                'FontSize', 9.5, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.2, 0.2, 0.2]);
        end

        m_u = sections{sec, 2}; s_u = sections{sec, 3};
        m_l = sections{sec, 4}; s_l = sections{sec, 5};
        m_e = sections{sec, 6}; s_e = sections{sec, 7};

        for a = 1:n_algs
            y_cursor = y_cursor - row_h;

            if mod(a, 2) == 0
                row_bg = [0.96, 0.96, 0.98];
            else
                row_bg = [1, 1, 1];
            end
            rectangle('Position', [col_x(1) - 10, y_cursor - row_h/2, sum(col_widths_px) + 20, row_h], ...
                'FaceColor', row_bg, 'EdgeColor', [0.88, 0.88, 0.88], 'Curvature', 0);

            text(ax, col_x(1) + 8, y_cursor, labels{a}, ...
                'FontSize', 9.5, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                'Color', colors(a, :));

            text(ax, col_x(2) + col_widths_px(2)/2, y_cursor, sprintf('%.2f \\pm %.2f', m_u(a), s_u(a)), ...
                'FontSize', 9.5, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.15, 0.15, 0.15]);

            text(ax, col_x(3) + col_widths_px(3)/2, y_cursor, sprintf('%.4f \\pm %.4f', m_l(a), s_l(a)), ...
                'FontSize', 9.5, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.15, 0.15, 0.15]);

            text(ax, col_x(4) + col_widths_px(4)/2, y_cursor, sprintf('%.2f \\pm %.2f', m_e(a), s_e(a)), ...
                'FontSize', 9.5, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.15, 0.15, 0.15]);
        end

        y_cursor = y_cursor - 10;
    end

    saveas(fig, fullfile(output_dir, 'extreme_solution_comparison.fig'));
    saveas(fig, fullfile(output_dir, 'extreme_solution_comparison.png'));
    close(fig);
end

function drawEfficiencyRadar(labels, norm_ku, norm_kl, norm_ke, energy_efficiency, composite_score, colors, output_dir, scene)
    n_algs = length(labels);

    max_eff = max(energy_efficiency);
    min_eff = min(energy_efficiency);
    range_eff = max_eff - min_eff;
    if range_eff < 1e-10, range_eff = 1; end
    norm_eff = (energy_efficiency - min_eff) / range_eff;

    fig = figure('Position', [100, 100, 700, 600], 'Color', 'w');
    hold on;

    n_axes = 4;
    angles = linspace(0, 2*pi, n_axes + 1);
    angles = angles(1:end-1);

    axis_labels = {'Utility', 'Latency', 'Energy', 'U/E'};

    for ring = [0.2, 0.4, 0.6, 0.8, 1.0]
        ring_x = ring * cos(angles);
        ring_y = ring * sin(angles);
        ring_x = [ring_x, ring_x(1)];
        ring_y = [ring_y, ring_y(1)];
        plot(ring_x, ring_y, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);
    end

    for k = 1:n_axes
        plot([0, cos(angles(k))], [0, sin(angles(k))], '-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.8);
    end

    for a = 1:n_algs
        vals = [norm_ku(a), norm_kl(a), norm_ke(a), norm_eff(a)];
        x_pts = vals .* cos(angles);
        y_pts = vals .* sin(angles);
        x_pts = [x_pts, x_pts(1)];
        y_pts = [y_pts, y_pts(1)];

        if a == 1
            lw = 2.5;
            fa = 0.15;
        else
            lw = 1.2;
            fa = 0.05;
        end

        fill(x_pts, y_pts, colors(a, :), 'FaceAlpha', fa, 'EdgeColor', colors(a, :), 'LineWidth', lw, ...
            'DisplayName', sprintf('%s (D=%.3f)', labels{a}, composite_score(a)));
    end

    for k = 1:n_axes
        label_r = 1.18;
        text(label_r * cos(angles(k)), label_r * sin(angles(k)), axis_labels{k}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    axis equal;
    axis off;
    title(sprintf('Efficiency-Enhanced Radar (%s Scene)', scene), ...
        'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 14);
    legend('Location', 'southoutside', 'FontName', 'Times New Roman', 'FontSize', 9, 'NumColumns', 2);

    saveas(fig, fullfile(output_dir, 'knee_point_efficiency_radar.fig'));
    saveas(fig, fullfile(output_dir, 'knee_point_efficiency_radar.png'));
    close(fig);
end
