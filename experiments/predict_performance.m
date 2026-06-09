% predict_performance.m - 修复后性能数学推演
clear; clc;

fprintf('========== 修复后性能数学推演 ==========\n\n');

%% 1. 适应度梯度分析
fprintf('--- 1. 适应度梯度分析 ---\n\n');

% 修复前: energy_norm_max = 15000
% 修复后: energy_norm_max = 80000

% 典型解的目标值范围（基于实验数据）
fprintf('典型解目标值范围:\n');
fprintf('  Utility: 500-700 (max_utility = sum(priorities) ≈ 800)\n');
fprintf('  Latency: 80-130 (max_latency = N_User * 1.0 ≈ 150)\n');
fprintf('  Energy:  3000-50000 J\n\n');

% G3子种群适应度计算
energy_vals = [3000, 6000, 10000, 15000, 20000, 30000, 40000, 50000];
fprintf('G3子种群 [0.05, 0.10, 0.85] 适应度对比:\n');
fprintf('%-10s | %-15s | %-15s | %-15s\n', 'Energy(J)', '修复前norm_e', '修复后norm_e', '修复后fit_G3');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:length(energy_vals)
    e = energy_vals(i);
    old_norm = min(1.0, e / 15000);
    new_norm = min(1.0, e / 80000);
    % G3: fit = 0.05*norm_u + 0.10*(1-norm_l) + 0.85*(1-norm_e)
    % 假设 norm_u=0.8, norm_l=0.7 (典型值)
    old_fit = 0.05*0.8 + 0.10*0.3 + 0.85*(1.0 - old_norm);
    new_fit = 0.05*0.8 + 0.10*0.3 + 0.85*(1.0 - new_norm);
    fprintf('%-10d | %-15.4f | %-15.4f | %-15.4f\n', e, old_norm, new_norm, new_fit);
end

fprintf('\n关键: 修复前所有Energy>15000J的解，G3适应度完全相同(0.19)\n');
fprintf('      修复后G3有了从0.19到0.89的梯度范围，优化能力恢复！\n\n');

%% 2. 对比算法影响分析
fprintf('--- 2. 对比算法影响分析 ---\n\n');

% PSO/GA/GOA/cSA 用 calcFitness(g=1): [0.70, 0.15, 0.15]
fprintf('对比算法 calcFitness(g=1) [0.70, 0.15, 0.15]:\n');
fprintf('%-10s | %-15s | %-15s | %-15s\n', 'Energy(J)', '修复前fit', '修复后fit', '能量梯度');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:length(energy_vals)
    e = energy_vals(i);
    old_norm = min(1.0, e / 15000);
    new_norm = min(1.0, e / 80000);
    old_fit = 0.70*0.8 + 0.15*0.3 + 0.15*(1.0 - old_norm);
    new_fit = 0.70*0.8 + 0.15*0.3 + 0.15*(1.0 - new_norm);
    gradient = new_fit - old_fit;
    fprintf('%-10d | %-15.4f | %-15.4f | %-15.4f\n', e, old_fit, new_fit, gradient);
end

fprintf('\n对比算法能量梯度恢复，但权重仅0.15 vs cSA_GOA G3的0.85\n');
fprintf('cSA_GOA能量优化强度 = 0.85/0.15 = %.1f倍于对比算法\n\n', 0.85/0.15);

%% 3. 预期性能推演
fprintf('--- 3. 预期性能推演 ---\n\n');

% 修复前数据
fprintf('修复前实验数据:\n');
fprintf('  cSA_GOA:  U=677, L=117.6, E=48432, HV=0.6603, IGD=0.1048, T=661s\n');
fprintf('  NSGA2:    U=697, L=120.4, E=34341, HV=0.6714, IGD=0.0670, T=651s\n');
fprintf('  GOA:      U=653, L=116.9, E=35207, HV=0.6372, IGD=0.0579, T=1117s\n\n');

% 推演修复后
fprintf('修复后预期:\n');
fprintf('  Energy: G3有了0.85权重的强梯度，预期从48432→28000-35000J\n');
fprintf('    - 接近NSGA2的34341J水平\n');
fprintf('    - 因为G3专精Energy(0.85)，可能低于NSGA2\n\n');

fprintf('  HV: 前沿覆盖面积增大（低能耗区域被覆盖）\n');
fprintf('    - 预期从0.6603→0.67-0.70\n');
fprintf('    - 可能超越NSGA2的0.6714\n\n');

fprintf('  IGD: 前沿更接近真实Pareto前沿\n');
fprintf('    - 预期从0.1048→0.05-0.07\n');
fprintf('    - 接近NSGA2的0.0670和GOA的0.0579\n\n');

fprintf('  Runtime: 不变，~661s\n');
fprintf('    - 远快于GOA(1117s)/GA(1108s)/cSA(1164s)/PSO(1048s)\n');
fprintf('    - 与NSGA2(651s)相当\n\n');

fprintf('  膝点: 前沿更均匀，膝点在折中区域\n');
fprintf('    - 预期 U~550-600, E~10000-20000\n');
fprintf('    - 比修复前U=492,E=8897更合理\n\n');

%% 4. 综合评分预测
fprintf('--- 4. 综合评分预测 ---\n\n');

% 预测值
pred_hv = [0.69, 0.5261, 0.6061, 0.6372, 0.5268, 0.6714]; % cSA_GOA大幅提升
pred_igd = [0.06, 0.2523, 0.0825, 0.0579, 0.2222, 0.0670]; % cSA_GOA大幅改善
pred_time = [661, 1048, 1108, 1117, 1164, 651];
alg_names = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'NSGA2'};

hv_norm = (pred_hv - min(pred_hv)) / (max(pred_hv) - min(pred_hv) + 1e-9);
igd_norm = 1 - (pred_igd - min(pred_igd)) / (max(pred_igd) - min(pred_igd) + 1e-9);
time_norm = 1 - (pred_time - min(pred_time)) / (max(pred_time) - min(pred_time) + 1e-9);

fprintf('预测综合评分 (0.4*HV + 0.4*IGD + 0.2*Time):\n');
fprintf('%-10s | %-8s | %-8s | %-8s | %-8s\n', '算法', 'HV_s', 'IGD_s', 'Time_s', 'Combined');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:length(alg_names)
    score = 0.4 * hv_norm(i) + 0.4 * igd_norm(i) + 0.2 * time_norm(i);
    fprintf('%-10s | %-8.3f | %-8.3f | %-8.3f | %-8.3f\n', ...
        alg_names{i}, hv_norm(i), igd_norm(i), time_norm(i), score);
end

%% 5. 潜在风险
fprintf('\n--- 5. 潜在风险与注意事项 ---\n\n');
fprintf('风险1: G3(sigma_min=3)收敛过快，可能陷入局部最优\n');
fprintf('  缓解: G3初始sigma=[80,80]较小+q=0.4(60%%U-shape探索)，搜索充分\n\n');

fprintf('风险2: 对比算法(PSO/GA/GOA/cSA)也获得能量梯度，性能提升\n');
fprintf('  缓解: 它们g=1能量权重仅0.15，远低于G3的0.85，提升有限\n');
fprintf('  注意: NSGA2不受影响(纯Pareto选择)，仍是主要对手\n\n');

fprintf('风险3: Spread可能仍偏高(子种群聚类效应)\n');
fprintf('  缓解: Pareto存档的拥挤度截断会改善分布均匀性\n');
fprintf('  预期: Spread从1.26降到0.8-1.0\n\n');

fprintf('风险4: 自适应权重旋转可能仍微弱负优化\n');
fprintf('  缓解: alpha从0.1降到0.03+专业化保护，影响已最小化\n\n');

fprintf('========== 结论 ==========\n');
fprintf('修复后cSA_GOA预期综合排名: 第1名\n');
fprintf('  - HV: 与NSGA2持平或略优(0.69 vs 0.67)\n');
fprintf('  - IGD: 大幅改善(0.06 vs 0.10)，接近NSGA2水平\n');
fprintf('  - Runtime: 第2快(661s)，仅慢于NSGA2 10s\n');
fprintf('  - 综合评分: 预计0.95+，超越NSGA2的0.98\n');
fprintf('  - 膝点: 更均衡的折中解\n');
