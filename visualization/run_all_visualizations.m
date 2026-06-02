function run_all_visualizations()
% run_all_visualizations - 对所有场景（Small/Medium/Large）运行完整可视化
% 自动识别 experiments 文件夹中的 mat 文件，按场景分类并生成可视化
fprintf('========== Running All Visualizations ==========\n\n');

scenes = {'Small', 'Medium', 'Large'};

for s = 1:length(scenes)
    scene = scenes{s};
    fprintf('========== Processing Scene: %s ==========\n', scene);

    comparison_pattern = fullfile('..', 'experiments', ['comparison_results_para_*' scene '_*.mat']);
    comp_files = dir(comparison_pattern);
    if ~isempty(comp_files)
        [~, idx] = sort([comp_files.datenum], 'descend');
        latest_comp = fullfile(comp_files(idx(1)).folder, comp_files(idx(1)).name);
        fprintf('  Latest comparison file: %s\n', comp_files(idx(1)).name);
        try
            plot_comparison_all(latest_comp);
            fprintf('   -> Comparison charts done\n');
        catch ME
            fprintf('   -> Error: %s\n', ME.message);
        end

        try
            plot_knee_comparison(latest_comp);
            fprintf('   -> Knee point comparison done\n\n');
        catch ME
            fprintf('   -> Knee point error: %s\n\n', ME.message);
        end
    else
        fprintf('  No comparison results found for %s scene\n\n', scene);
    end

    ablation_pattern = fullfile('..', 'experiments', ['ablation_results_para_*' scene '_*.mat']);
    abl_files = dir(ablation_pattern);
    if ~isempty(abl_files)
        [~, idx] = sort([abl_files.datenum], 'descend');
        latest_abl = fullfile(abl_files(idx(1)).folder, abl_files(idx(1)).name);
        fprintf('  Latest ablation file: %s\n', abl_files(idx(1)).name);
        try
            plot_ablation_all(latest_abl);
            fprintf('   -> Ablation charts done\n\n');
        catch ME
            fprintf('   -> Error: %s\n\n', ME.message);
        end
    else
        fprintf('  No ablation results found for %s scene\n\n', scene);
    end
end

fprintf('========== All Visualizations Complete ==========\n');
fprintf('Please check the following directories for your figures:\n');
fprintf('  - figures/comparison/Small\n');
fprintf('  - figures/comparison/Medium\n');
fprintf('  - figures/comparison/Large\n');
fprintf('  - figures/ablation/Small\n');
fprintf('  - figures/ablation/Medium\n');
fprintf('  - figures/ablation/Large\n');
end
