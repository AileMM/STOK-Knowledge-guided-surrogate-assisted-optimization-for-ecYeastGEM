function [Obj, Con] = bio10(var, M, ~)
    
    [ps, D] = size(var);  % 种群个体数，变量维度
    Obj = zeros(ps, M);

    data = load('sourceModel_C_turanose_N_L_glutamate.mat');
    kcats_inf_struct = load('kcats_inf.mat');
    
    kcats_inf_cell = kcats_inf_struct.kcats_inf;
    oroginal_model = data.ecModel_C_turanose_N_L_glutamate; 
    ATP_row = 324;  % ATP在S矩阵中的行编号（S_0434）作为反应物参与
    changeCobraSolver('gurobi', 'LP', 0);
    for i = 1:ps
        model = oroginal_model; % 复制模型
%         origin_sol = optimizeCbModel(model, 'max');
        for j = 1:D
            rc = kcats_inf_cell{j};
            row_idx = rc{1};
            col_idx = rc{2};
            model.S(row_idx, col_idx) = -var(i, j);  
        end
        sol = optimizeCbModel(model, 'max'); 
        ref_val = 0.3953;  
        F1 = abs(ref_val - sol.f); 
        
        % F2: ATP消耗（第324行中所有反应为ATP反应物的通量之和）
        ATP_flux = 0;
        for rxn = 1:length(sol.v)
            if model.S(ATP_row, rxn) < 0
                ATP_flux = ATP_flux + abs(sol.v(rxn));
            end
        end
        F2 = ATP_flux;
        Obj(i, :) = [F1, F2];
    end
end
