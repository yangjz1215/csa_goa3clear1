% calcCapturability.m - 捕获能力计算函数（按设计方案修改）
% 公式：Capturability_g^t = 1/(R·t₂)
% 其中：R = (M·Vel²)/L（塘鹅动能系数），t₂ = 1 + It/T_max（时间衰减项）
function capturability = calcCapturability(subpop, iter, FES_max, g)
    % 捕获能力计算（子种群差异化）
    if g == 3
        % G3（能耗优化）捕获能力更强，更早收敛
        capturability = 0.3 + 0.6 * (iter / FES_max);  % 0.3~0.9
    else
        % G1和G2捕获能力
        capturability = 0.1 + 0.8 * (iter / FES_max);  % 0.1~0.9
    end
    
    % 确保在[0,1]范围内
    capturability = max(0, min(1, capturability));
end