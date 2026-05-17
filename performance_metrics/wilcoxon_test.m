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
