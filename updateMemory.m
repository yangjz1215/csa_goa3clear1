function mem = updateMemory(old_mem, new_candidates, User, priorities, E_remaining, E_max, ...
    k_move, g, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params)
%UPDATEMEMORY 更新记忆矩阵：按加权适应度排序选优，同时对子目标最优解做配额保护
%   配额保护确保效用/时延/能耗各自的极端解不被淘汰，维持Pareto前沿覆盖广度

    K_old = size(old_mem, 1);
    K_new = size(new_candidates, 1);
    all_candidates = zeros(K_old + K_new, N_UAV, 2);
    all_candidates(1:K_old, :, :) = old_mem;
    all_candidates(K_old+1:end, :, :) = new_candidates;

    if nargin < 13 || isempty(capturability)
        capturability = 0.5;
    end

    if isfield(params, 'variant') && strcmp(params.variant, 'no_subpop')
        g_eval = 1;
    else
        g_eval = g;
    end

    n = K_old + K_new;
    fits = zeros(1, n);
    utils = zeros(1, n);
    lats = zeros(1, n);
    nrgs = zeros(1, n);

    % 一次循环同时获取加权适应度和子目标值，避免重复调用calcMEC_Objectives
    for i = 1:n
        candidate = squeeze(all_candidates(i, :, :));
        if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
            candidate = reshape(candidate, N_UAV, 2);
        end
        [fits(i), utils(i), lats(i), nrgs(i)] = calcFitness(candidate, User, priorities, ...
            E_remaining, E_max, k_move, g_eval, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params);
    end

    % 配额保护：每子目标保留m个最优解，不被加权适应度淘汰
    m_quota = 0;
    if isfield(params, 'mem_quota_m') && ~isempty(params.mem_quota_m) && params.mem_quota_m > 0 && K_old > 0
        m_quota = min(round(params.mem_quota_m), max(1, floor(K_old / 3)));
    end

    if m_quota <= 0
        [~, idx] = sort(fits, 'descend');
        mem = all_candidates(idx(1:K_old), :, :);
        return;
    end

    % 按当前子群对应的子目标排序，保护极端解
    switch g_eval
        case 1
            [~, ord_lead] = sort(utils, 'descend');
        case 2
            [~, ord_lead] = sort(lats, 'ascend');
        otherwise
            [~, ord_lead] = sort(nrgs, 'ascend');
    end

    prot = ord_lead(1:min(m_quota, n));
    rest = setdiff(1:n, prot, 'stable');
    need = K_old - numel(prot);
    if need <= 0
        mem = all_candidates(prot(1:K_old), :, :);
        return;
    end

    % 剩余位置按加权适应度填充
    [~, ord_fit] = sort(fits(rest), 'descend');
    take = min(need, numel(rest));
    idx_rest = rest(ord_fit(1:take));
    chosen = [prot(:)', idx_rest];
    mem = all_candidates(chosen, :, :);
end
