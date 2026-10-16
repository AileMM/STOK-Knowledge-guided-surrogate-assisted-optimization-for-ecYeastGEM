function checkMediaChange(model, C_source, N_source)
    % 检查碳源exchange上下界
    fprintf('Checking Carbon source: %s\n', C_source);
    [lbC, ubC] = getExchangeBounds(model, C_source);
    if isempty(lbC)
        fprintf('  Warning: Carbon source exchange reaction not found!\n');
    else
        fprintf('  Lower bound (lb): %g\n', lbC);
        fprintf('  Upper bound (ub): %g\n', ubC);
    end
    
    % 检查氮源exchange上下界
    fprintf('Checking Nitrogen source: %s\n', N_source);
    [lbN, ubN] = getExchangeBounds(model, N_source);
    if isempty(lbN)
        fprintf('  Warning: Nitrogen source exchange reaction not found!\n');
    else
        fprintf('  Lower bound (lb): %g\n', lbN);
        fprintf('  Upper bound (ub): %g\n', ubN);
    end
end

function [lb, ub] = getExchangeBounds(model, source)
    % 返回某种物质exchange和exchange (reversible)反应的lb和ub
    lb = [];
    ub = [];
    rxnName_rev = [source ' exchange (reversible)'];
    rxnName_fwd = [source ' exchange'];
    
    pos_rev = find(strcmp(model.rxnNames, rxnName_rev));
    pos_fwd = find(strcmp(model.rxnNames, rxnName_fwd));
    
    if ~isempty(pos_rev)
        lb = model.lb(pos_rev);
        ub = model.ub(pos_rev);
    elseif ~isempty(pos_fwd)
        lb = model.lb(pos_fwd);
        ub = model.ub(pos_fwd);
    end
end
