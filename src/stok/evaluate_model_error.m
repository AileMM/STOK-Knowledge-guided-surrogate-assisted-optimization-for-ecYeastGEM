function rmse = evaluate_model_error(func_surrogate, database,M)

    Out_dec = database(:,1:end-M); 
    Out_obj = database(:,end-M+1:end);
    % 随机取一小部分数据作为测试集
    idx = randperm(size(Out_dec,1), min(20, size(Out_dec,1)));
    X_test = Out_dec(idx,:);
    Y_true = Out_obj(idx,:);
    
    % 代理模型预测
     for i = 1:M
        [objs_pre(:, i),objs_std(:, i)] = func_surrogate.gpr{i}.func(X_test);
    end
    Y_pred = objs_pre;
    
    % 计算RMSE
    rmse = sqrt(mean((Y_true - Y_pred).^2, 'all'));
end
