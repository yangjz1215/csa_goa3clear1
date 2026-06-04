function [archive, is_updated] = updateParetoArchive3D(archive, new_pos, util, lat, nrg)
    is_updated = false;
    new_obj = [-util, lat, nrg];

    is_dominated = false;
    indices_to_remove = [];

    for i = 1:length(archive)
        old_obj = [-archive(i).Utility, archive(i).Latency, archive(i).Energy];

        if all(old_obj <= new_obj) && any(old_obj < new_obj)
            is_dominated = true;
            break;
        end

        if all(new_obj <= old_obj) && any(new_obj < old_obj)
            indices_to_remove = [indices_to_remove, i];
        end
    end

    if ~is_dominated
        archive(indices_to_remove) = [];
        archive(end+1) = struct('UAV_pos', new_pos, 'Utility', util, 'Latency', lat, 'Energy', nrg);
        is_updated = true;
    end

    MAX_ARCHIVE_SIZE = 500;
    % 性能优化：允许档案临时超限到1.2倍再清理，避免每次新插入都触发O(n^2)的全局支配清理
    % 原代码每次超限都调用 removeDominated(archive)，存档满时是25万次比较，性能瓶颈
    if length(archive) > MAX_ARCHIVE_SIZE * 1.2
        % 截断前先做一次全局支配清理，确保档案中不含被支配解
        archive = removeDominated(archive);

        if length(archive) <= MAX_ARCHIVE_SIZE
            return;
        end

        arch_U = [archive.Utility];
        arch_L = [archive.Latency];
        arch_E = [archive.Energy];

        [~, max_u_idx] = max(arch_U);
        [~, min_l_idx] = min(arch_L);
        [~, min_e_idx] = min(arch_E);
        extreme_indices = unique([max_u_idx, min_l_idx, min_e_idx]);

        all_indices = 1:length(archive);
        other_indices = setdiff(all_indices, extreme_indices);

        num_to_keep = MAX_ARCHIVE_SIZE - length(extreme_indices);
        if num_to_keep <= 0
            archive = archive(extreme_indices);
        elseif isempty(other_indices)
            archive = archive(extreme_indices);
        elseif length(other_indices) <= num_to_keep
            final_keep_indices = [extreme_indices, other_indices];
            archive = archive(final_keep_indices);
        else
            keep_others = selectArchiveIndicesByCrowding3D(archive, other_indices, num_to_keep);
            final_keep_indices = [extreme_indices, keep_others];
            archive = archive(final_keep_indices);
        end
    end
end

function keep_idx = selectArchiveIndicesByCrowding3D(archive, other_indices, num_to_keep)
% Crowding-distance truncation (Deb et al.) on 3 objectives (all minimization):
%   f1 = -Utility, f2 = Latency, f3 = Energy
    n = numel(other_indices);
    if n <= num_to_keep
        keep_idx = other_indices(:)';
        return;
    end

    F = zeros(n, 3);
    for ii = 1:n
        a = archive(other_indices(ii));
        F(ii, :) = [-a.Utility, a.Latency, a.Energy];
    end

    cd = zeros(n, 1);
    for m = 1:3
        [sorted_f, ord] = sort(F(:, m));
        cd(ord(1)) = inf;
        cd(ord(end)) = inf;
        span = max(sorted_f(end) - sorted_f(1), 1e-12);
        for j = 2:n - 1
            idx_j = ord(j);
            cd(idx_j) = cd(idx_j) + (sorted_f(j + 1) - sorted_f(j - 1)) / span;
        end
    end

    [~, rank] = sort(cd, 'descend');
    pick = rank(1:num_to_keep);
    keep_idx = other_indices(pick);
end

function archive = removeDominated(archive)
% 全局支配清理：移除档案中所有被支配的解
    n = length(archive);
    dominated = false(n, 1);
    for i = 1:n
        if dominated(i), continue; end
        oi_i = [-archive(i).Utility, archive(i).Latency, archive(i).Energy];
        for j = i+1:n
            if dominated(j), continue; end
            oi_j = [-archive(j).Utility, archive(j).Latency, archive(j).Energy];
            if all(oi_i <= oi_j) && any(oi_i < oi_j)
                dominated(j) = true;
            elseif all(oi_j <= oi_i) && any(oi_j < oi_i)
                dominated(i) = true;
                break;
            end
        end
    end
    archive = archive(~dominated);
end
