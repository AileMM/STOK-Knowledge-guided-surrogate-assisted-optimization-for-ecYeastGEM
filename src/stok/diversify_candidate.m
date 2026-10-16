function new_candidate = diversify_candidate(solution_candidate, archive_best_dec, lb, ub)
    epsilon = 5e-3;
    no_trials = 50;
    scales = linspace(0.1, 1, no_trials);
    c = 1;
    dim = size(solution_candidate, 2);
    new_candidate = solution_candidate;

    for i = 1:size(solution_candidate, 1)
        while min(max(abs(repmat(new_candidate(i,:), size(archive_best_dec,1), 1) - ...
                archive_best_dec), [], 2)) < epsilon
            scale = scales(mod(c-1,no_trials)+1);
            perturbation = scale * (ub - lb) .* (rand(1, dim) - 0.5);
            new_candidate(i,:) = new_candidate(i,:) + perturbation;
            new_candidate(i,:) = min(max(new_candidate(i,:), lb), ub);
            c = c + 1;
        end
    end
end
