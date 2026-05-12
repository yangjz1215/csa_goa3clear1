function p = norm_cdf_standard(x)
%NORM_CDF_STANDARD Standard normal CDF without Statistics Toolbox.
    p = 0.5 * (1 + erf(x ./ sqrt(2)));
end
