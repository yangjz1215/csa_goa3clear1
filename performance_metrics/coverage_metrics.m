function coverage_metrics = calculate_coverage_metrics(UAV_pos, User, priorities, params, RRH)
    if nargin < 6
        RRH = [];
    end

    N_User = size(User, 1);
    cover_radius = params.cover_radius;
    RRH_radius = params.RRH_radius;

    covered_by_uav = false(N_User, 1);
    covered_by_rrh = false(N_User, 1);
    covered_count = 0;
    high_priority_covered = 0;
    total_priority_sum = 0;
    covered_priority_sum = 0;

    for i = 1:N_User
        dists_uav = sqrt(sum((UAV_pos - User(i, :)).^2, 2));
        covered_by_uav(i) = any(dists_uav <= cover_radius);

        if ~isempty(RRH)
            dists_rrh = sqrt(sum((RRH - User(i, :)).^2, 2));
            covered_by_rrh(i) = any(dists_rrh <= RRH_radius);
        end

        if covered_by_uav(i) || covered_by_rrh(i)
            covered_count = covered_count + 1;
            covered_priority_sum = covered_priority_sum + priorities(i);

            if priorities(i) >= 3
                high_priority_covered = high_priority_covered + 1;
            end
        end

        total_priority_sum = total_priority_sum + priorities(i);
    end

    coverage_metrics.total_coverage_ratio = covered_count / N_User;
    coverage_metrics.high_priority_coverage_ratio = high_priority_covered / sum(priorities >= 3);
    coverage_metrics.priority_coverage_ratio = covered_priority_sum / total_priority_sum;
    coverage_metrics.covered_user_count = covered_count;
    coverage_metrics.high_priority_covered_count = high_priority_covered;
end
