function [Obj, Con] = bio1(var, M, ~)
    % bio1 多目标优化函数
    % var: N x D 的个体矩阵，每行一个个体，对应所有 kcat
    % M: 目标数，此处为 2
    % 输出:
    %   Obj: N x M 矩阵，每行为一个个体的目标值 [误差, ATP 消耗]
    %   Con: 约束，此处为全零、
    % sol.f,最优目标值（即最大化的 biomass 通量）
    % sol.x 优化得到的所有反应的通量向量（flux vector）
    % sol.v最大速率
    % 但是注意画图的时候要考虑原模型的帕累托点
    
    [ps, D] = size(var);  % 种群个体数，变量维度
    Obj = zeros(ps, M);
    data = load('sourceModel_C__R__lactate_N_L_aspartate.mat');
    kcats_inf_struct = load('kcats_inf.mat');
    kcats_inf_cell = kcats_inf_struct.kcats_inf;
%     oroginal_model = data.ecModel_C__R__lactate_N_L_aspartate; 
    ATP_row = 324;  % ATP在S矩阵中的行编号（S_0434）
    %% 如果报错的话，需要重新添加cobratool的路径
    changeCobraSolver('gurobi', 'LP', 0);
    for i = 1:ps
        model_new = data.ecModel_C__R__lactate_N_L_aspartate; % 复制模型
%         origin_sol = optimizeCbModel(model_new, 'max');
        for j = 1:D
            rc = kcats_inf_cell{j};
            row_idx = rc{1};
            col_idx = rc{2};
            model_new.S(row_idx, col_idx) = -var(i, j);  %这里负号是反应物
        end
        sol = optimizeCbModel(model_new, 'max'); % 最大化生物量（biomass）反应的通量，最大生长率
        ref_val = 0.4123;  
        F1 = abs(ref_val - sol.f); 
        
        % F2: ATP消耗（第324行中所有反应为ATP反应物的通量之和）
        ATP_flux = 0;
        for rxn = 1:length(sol.v)
            if model_new.S(ATP_row, rxn) < 0
                ATP_flux = ATP_flux + abs(sol.v(rxn));
            end
        end
        F2 = ATP_flux;
        Obj(i, :) = [F1, F2];
    end
end
