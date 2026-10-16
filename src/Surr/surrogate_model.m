function [func_surrogate,model_type] = surrogate_model(X,y,surrogate_name,M)
no_samples = size(X,M);
dim = size(X,2);

%% 多个代理模型的构建
switch(surrogate_name)
    case 'prs'
        glmModel = fitglm(X,y,'purequadratic');
        func_surrogate = @(x)predict(glmModel,x);
        model_type = 'regression';
        
    case 'gpr'  %fitrgp 只支持单目标回归

%         gprModel = fitrgp(X, y, ...
%                 'Basis', 'constant', ...  % 使用常数基函数
%                 'KernelFunction', 'squaredexponential', ...  % 使用平方指数核
%                 'FitMethod', 'exact', ...
%                 'PredictMethod', 'exact', ...
%                 'Standardize', 1, ...
%                 'Sigma', 1e-6);  % 增加噪声项，避免矩阵奇异
% 
%         func_surrogate = @(x)predict(gprModel,x);
%         model_type = 'regression';
        
                   gprModel = fitrgp(X, y, ...
                'Basis', 'constant', ...
                'KernelFunction', 'squaredexponential', ...
                'FitMethod', 'exact', ...
                'PredictMethod', 'exact', ...
                'Standardize', 1, ...
                'Sigma', 1e-6);  % 增加噪声项，避免矩阵奇异

            func_surrogate = struct( ...
                'func', @(x)predict(gprModel, x), ...
                'model', gprModel, ...
                'type', 'gpr');

            model_type = 'regression';

        
    case 'rbf'
%         R = zeros(no_samples,no_samples);
%         for i=1:no_samples
%             for j=1:no_samples
%                 R(i,j)=sqrt(sum((X(i,:)-X(j,:)).^2));
%             end
%         end
%         spr_estimate = max(max(R))/(dim*no_samples)^(1/dim);
%         net = newrbe(X',y',spr_estimate);
%         func_surrogate = @(x)sim(net,x');
%         model_type = 'regression';
        % flag1核函数
        flag = 'cubic';
        [lambda1, gamma1] = RBF(X,y,flag,M);
        func_surrogate = @(x)RBF_eval(x,X,lambda1,gamma1,flag,M);
        
        
    case 'kriging'
        % 创建 Kriging/Gaussian Process 回归模型
        gprMdl = fitrgp(X, y, 'KernelFunction', 'ardsquaredexponential', ...
                        'BasisFunction', 'constant', 'Sigma', 1e-3, ...
                        'FitMethod', 'exact', 'PredictMethod', 'exact');

        % 代理模型预测函数
        func_surrogate = @(x) predict(gprMdl, x);
        % 模型类型
        model_type = 'regression';

end
