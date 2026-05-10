function hv = hypervolume(objectives, reference_point)
    if nargin < 2
        reference_point = max(objectives, [], 1) * 1.1;
    end

    n_points = size(objectives, 1);
    n_objectives = size(objectives, 2);

    if n_points == 0
        hv = 0;
        return;
    end

    obj_norm = normalize(objectives, reference_point);

    if n_objectives == 2
        hv = hv_2d(obj_norm);
    elseif n_objectives == 3
        hv = hv_3d(obj_norm);
    else
        hv = hv_nd(obj_norm, n_objectives);
    end
end

function obj_norm = normalize(objectives, reference_point)
    obj_norm = objectives;
    for m = 1:size(objectives, 2)
        obj_norm(:, m) = reference_point(m) - objectives(:, m);
    end
end

function hv = hv_2d(objectives)
    [~, idx] = sort(objectives(:, 1), 'descend');
    objectives = objectives(idx, :);

    hv = 0;
    prev_y = 0;

    for i = 1:size(objectives, 1)
        width = objectives(i, 1);
        height = objectives(i, 2) - prev_y;
        hv = hv + width * height;
        prev_y = objectives(i, 2);
    end
end

function hv = hv_3d(objectives)
    n = size(objectives, 1);
    if n == 0
        hv = 0;
        return;
    end

    [~, idx] = sort(objectives(:, 3), 'descend');
    objectives = objectives(idx, :);

    hv = 0;
    prev_z = 0;

    for i = 1:n
        z = objectives(i, 3);
        dz = z - prev_z;

        if dz > 0
            slice = objectives(i:end, 1:2);
            area = hv_2d(slice);
            hv = hv + area * dz;
        end

        prev_z = z;
    end
end

function hv = hv_nd(objectives, n_objectives)
    [~, idx] = sortrows(objectives, -1);
    objectives = objectives(idx, :);

    hv = 0;

    for i = 1:size(objectives, 1)
        if i == 1
            volume = prod(objectives(i, :));
            hv = hv + volume;
        else
            for j = 1:i-1
                diff = objectives(i, :) - objectives(j, :);
                if all(diff <= 0)
                    volume = prod(objectives(i, :));
                    hv = hv + volume;
                    break;
                end
            end
        end
    end
end
