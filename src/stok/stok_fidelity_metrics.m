function metrics = stok_fidelity_metrics(fidelity, rows)
% Per-objective metrics for paired, pre-evaluation predictions.
% rows optionally selects iterations, phases, or FE windows.
    if nargin < 2
        rows = true(size(fidelity.exact,1),1);
    end
    observed = fidelity.exact(rows,:);
    predicted = fidelity.predicted(rows,:);
    assert(isequal(size(observed),size(predicted)), ...
        'STOK:PairShape', 'Exact and predicted matrices must have identical sizes.');
    M = size(observed,2);
    N = zeros(M,1);
    Excluded = zeros(M,1);
    MAE = nan(M,1);
    RMSE = nan(M,1);
    R2 = nan(M,1);
    for j = 1:M
        valid = isfinite(observed(:,j)) & isfinite(predicted(:,j));
        N(j) = sum(valid);
        Excluded(j) = sum(~valid);
        if N(j) == 0
            continue;
        end
        y = observed(valid,j);
        residual = predicted(valid,j) - y;
        MAE(j) = mean(abs(residual));
        RMSE(j) = sqrt(mean(residual.^2));
        SST = sum((y - mean(y)).^2);
        if N(j) >= 2 && SST > 0
            R2(j) = 1 - sum(residual.^2)/SST;
        end
    end
    Objective = (1:M)';
    metrics = table(Objective,N,Excluded,MAE,RMSE,R2);
end
