function objs_acq_offspring = sur_Evaluation(func_surrogate,population_offspring,acquisition,M)

    switch(acquisition)
        case 'plain'
            objs_acq_offspring = func_surrogate(population_offspring);

        case 'lcb'
            w = 2;
            objs_pre = zeros(size(population_offspring, 1), M);
            objs_std = zeros(size(population_offspring, 1), M);
            if isfield(func_surrogate, 'rbf')
                objs_pre = func_surrogate.rbf(population_offspring);
                rbf_train_errors = calculate_rbf_training_errors(func_surrogate.rbf, population_offspring);
                objs_std = rbf_train_errors;  % 将训练误差扩展为标准差      
            elseif isfield(func_surrogate, 'gpr')
                for i = 1:M
                    [objs_pre(:, i),objs_std(:, i)] = func_surrogate.gpr{i}.func(population_offspring);
                end
            else
                for i = 1:M
                    [objs_pre(:, i),objs_std(:, i)] = func_surrogate.kri{i}(population_offspring);
                end
            
            end
%             objs_acq_offspring = objs_pre;
            objs_acq_offspring = objs_pre-w*objs_std;

        case 'lcb-decay'
            w = 2-2*(1/(1+exp(5-20*FEsUsed/500)));
            [objs_pre,objs_std] = func_surrogate(population_offspring);
            objs_acq_offspring = objs_pre-w*objs_std;

    end

    function rbf_train_errors = calculate_rbf_training_errors(rbf_model, population)
            y_pred = rbf_model(population); 
            [N, M] = size(y_pred);
            y_actual = rand(N, M);  
            errors = y_pred - y_actual;
            % 估算标准差（例如，使用残差的标准差作为不确定性）
%             rbf_train_errors = zscore(errors, 0, 1);  
            rbf_train_errors = 0.05 * ones(N, M); % 暂时设置较小常量值
    end
end