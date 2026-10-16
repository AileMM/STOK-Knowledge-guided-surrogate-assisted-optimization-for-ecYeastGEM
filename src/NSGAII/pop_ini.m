function [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,func_surrogate,acquisition,M)
population = database(:,1:end-M);
objs = database(:,end-M+1:end);

switch (ini_strategy)
%     case 'new'
%         population = population(end-popsize+1:end,:);
%         objs_acq = objs(end-popsize+1:end);
%     case 'elite'
%         [~,idx] = sort(objs);
%         population = population(idx(1:popsize),:);
%         objs_acq = objs(idx(1:popsize));
%     case 'random'
%         population = lhsdesign_modified(popsize,lb,ub);
%         switch (acquisition)
%             case 'plain'
%                 objs_acq = func_surrogate(population);
%             case 'lcb'
%                 w = 2;
%                 [objs_pre,objs_std] = func_surrogate(population);
%                 objs_acq = objs_pre-w*objs_std;
%             case 'lcb-decay'
%                 w = 2-2*(1/(1+exp(5-20*size(database,1)/500)));
%                 [objs_pre,objs_std] = func_surrogate(population);
%                 objs_acq = objs_pre-w*objs_std;
%         end
    case 'elite+random'
%         [~,idx] = sort(objs);
%         population = [population(idx(1:floor(popsize/2)),:);...
%             lhsdesign_modified(popsize-floor(popsize/2),lb,ub)];
        switch (acquisition)
%             case 'plain'
%                 objs_acq = func_surrogate(population);
                
            case 'lcb'       
                w = 2;
                [rank,~,~] = NSGASort(objs);
                idx = rank(1:popsize);
                population = population(idx,:);
                objs = objs(idx,:);
                objs_pre = zeros(size(population, 1), M, 3);  %  3 代表三种代理模型
                objs_std = zeros(size(population, 1), M, 3);
                for i = 1:M
                    % GPR 代理模型预测
                    [objs_pre(:, i, 1), objs_std(:, i, 1)] = func_surrogate.gpr{i}.func(population);
                    % kriging 代理模型预测
                    [objs_pre(:, i, 2), objs_std(:, i, 2)] = func_surrogate.kri{i}(population);
                end
                % RBF 代理模型预测（一次性返回所有目标的预测值）
                objs_pre(:, :, 3) = func_surrogate.rbf(population);
                rbf_train_errors = calculate_rbf_training_errors(func_surrogate.rbf, population);
                objs_std(:, :, 3) = rbf_train_errors;  % 将训练误差扩展为标准差

                % 计算 LCB 采集函数
                objs_acq = objs_pre - w * objs_std;  
                         
%             case 'lcb-decay'
%                 w = 2-2*(1/(1+exp(5-20*size(database,1)/500)));
%                 [objs_pre,objs_std] = func_surrogate(population);
%                 objs_acq = objs_pre-w*objs_std;
        end
%     case 'best+randperm'
%         [~,idx] = min(objs);
%         seq = randperm(length(objs));
%         population = [population(idx,:);population(seq(1:popsize-1),:)];
%         objs_acq = [objs(idx);objs(seq(1:popsize-1))];
end
    % 通过残差来估算不确定性：可以基于 RBF 模型的残差（训练误差） 来估算每个预测值的不确定性。
    function rbf_train_errors = calculate_rbf_training_errors(rbf_model, population)
        y_pred = rbf_model(population); 
        y_actual = objs;  
        errors = y_pred - y_actual;
        % 估算标准差（例如，使用残差的标准差作为不确定性）
        rbf_train_errors = zscore(errors, 0, 1);  
    end
end


