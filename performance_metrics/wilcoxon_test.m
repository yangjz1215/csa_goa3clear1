function [p_value, h_stat] = wilcoxon_test(data1, data2)
%WILCOXON_TEST Two-sample Wilcoxon rank-sum (Mann-Whitney U), no Statistics Toolbox.
%   Uses normal approximation with tie correction (same ranks as ranksum).

    data1 = data1(:);
    data2 = data2(:);
    n1 = numel(data1);
    n2 = numel(data2);
    N = n1 + n2;

    if n1 < 1 || n2 < 1
        p_value = NaN;
        h_stat = NaN;
        return;
    end

    vals = [data1; data2];
    grp1 = [true(n1, 1); false(n2, 1)];
    [sv, ord] = sort(vals);
    r_sorted = average_ranks_sorted(sv);
    ranks = zeros(N, 1);
    ranks(ord) = r_sorted;

    R1 = sum(ranks(grp1));
    U1 = R1 - n1 * (n1 + 1) / 2;
    mu = n1 * n2 / 2;
    sum_r2 = sum(ranks .^ 2);
    denom = N * (N - 1);
    if denom <= 0
        p_value = 1;
        h_stat = 0;
        return;
    end

    var_U = (n1 * n2 / denom) * (sum_r2 - N * (N + 1) ^ 2 / 4);
    if var_U <= 0 || ~isfinite(var_U)
        p_value = 1;
        h_stat = 0;
        return;
    end

    z = (U1 - mu) / sqrt(var_U);
    p_value = 2 * (1 - norm_cdf_standard(abs(z)));
    p_value = max(0, min(1, p_value));
    h_stat = double(p_value < 0.05);
end

function [p_value, h_stat] = wilcoxon_signed_rank_test(data1, data2)
    diff = data1 - data2;
    diff = diff(diff ~= 0);
    n = length(diff);

    if n == 0
        p_value = 1;
        h_stat = 0;
        return;
    end

    [abs_diff, idx] = sort(abs(diff));
    ranks = tie_correction(abs_diff);

    signed_ranks = zeros(n, 1);
    for i = 1:n
        if diff(idx(i)) > 0
            signed_ranks(i) = ranks(i);
        else
            signed_ranks(i) = -ranks(i);
        end
    end

    W = sum(signed_ranks);
    sigma = sqrt(n * (n + 1) * (2 * n + 1) / 6);
    if sigma <= 0 || ~isfinite(sigma)
        p_value = 1;
    else
        z = abs(W) / sigma;
        p_value = 2 * (1 - norm_cdf_standard(z));
        p_value = max(0, min(1, p_value));
    end

    h_stat = double(p_value < 0.05);
end

function ranks = tie_correction(abs_diff)
    n = length(abs_diff);
    ranks = zeros(n, 1);
    i = 1;

    while i <= n
        j = i;
        while j < n && abs_diff(j + 1) == abs_diff(i)
            j = j + 1;
        end

        rank_sum = sum(i:j);
        avg_rank = rank_sum / (j - i + 1);
        ranks(i:j) = avg_rank;

        i = j + 1;
    end
end

function r = average_ranks_sorted(sv)
    n = numel(sv);
    r = zeros(n, 1);
    k = 1;
    while k <= n
        kt = k;
        while kt < n && sv(kt + 1) == sv(k)
            kt = kt + 1;
        end
        r(k:kt) = (k + kt) / 2;
        k = kt + 1;
    end
end

function [p_value, h_stat] = statistical_summary(data, alpha)
    if nargin < 2
        alpha = 0.05;
    end

    n = length(data);
    mean_val = mean(data);
    std_val = std(data);
    median_val = median(data);

    fprintf('样本数: %d\n', n);
    fprintf('均值: %.4f\n', mean_val);
    fprintf('标准差: %.4f\n', std_val);
    fprintf('中位数: %.4f\n', median_val);
    fprintf('最小值: %.4f\n', min(data));
    fprintf('最大值: %.4f\n', max(data));

    if nargout == 0
        fprintf('95%%置信区间: [%.4f, %.4f]\n', ...
            mean_val - 1.96 * std_val / sqrt(n), ...
            mean_val + 1.96 * std_val / sqrt(n));
    end

    p_value = NaN;
    h_stat = NaN;
end
