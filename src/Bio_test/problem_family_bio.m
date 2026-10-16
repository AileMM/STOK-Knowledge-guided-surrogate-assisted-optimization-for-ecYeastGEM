function problem= problem_family_bio(problem_id,M,d)
% 每个模型都要重新修改C，N源，fnc里面都是要重新更改模型参数的
 % 从 kcats_inf 中提取所有子 cell 的第3个值
data= load('kcats_inf.mat');
kcats_inf = data.kcats_inf;  % 取出真正的 cell 数组
all_kcat_values = cellfun(@(c) c{1,3}, kcats_inf, 'UniformOutput', false);  % 得到 cell 数组
all_kcat_values = cell2mat(all_kcat_values);  % 转成数值行向量
switch(problem_id)
    case 1 
        dim = d;
        problem.M = M;
        problem.name = 'bio1';
        problem.dims = dim; 
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        
        optimum = UniformPoint(500,M)/2;
        problem.fnc = @(x)bio1(x,problem.M,optimum);
        problem.F_best = optimum;
    case 2 
        dim = d;
        problem.M = M;
        problem.name = 'bio2';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        problem.fnc = @(x)bio2(x,problem.M,optimum);
        problem.F_best = optimum;
    case 3 
        dim = d;
        problem.M = M;
        problem.name = 'bio3';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        problem.fnc = @(x)bio3(x,problem.M,optimum);
        problem.F_best = optimum;
    case 4 
        dim = d;
        problem.M = M;
        problem.name = 'bio4';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        problem.fnc = @(x)bio4(x,problem.M,optimum);
        problem.F_best = optimum;
    case 5 
        dim = d;
        problem.M = M;
        problem.name = 'bio5';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        
        optimum = [0:1/(500-1):1;1:-1/(500-1):0]';
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,size(optimum,2));
        problem.fnc = @(x)bio5(x,problem.M,optimum);
        problem.F_best = optimum;
    case 6 
        dim = d;
        problem.M = M;
        problem.name = 'bio6';
        problem.dims = dim;
        problem.lb = 0.001 * all_kcat_values;
        problem.ub = 10 * all_kcat_values;
        optimum = [0:1/(500-1):1;1:-1/(500-1):0]';
        problem.fnc = @(x)bio6(x,problem.M,optimum);
        problem.F_best = optimum;
    case 7 
        dim = d;
        problem.M = M;
        problem.name = 'bio7';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = UniformPoint(500,M-1,'grid');
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        optimum            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
        rotation = eye(dim,dim);
        problem.fnc = @(x)bio7(x,problem.M,optimum);
        problem.F_best = optimum;
    case 8 
        dim = d;
        problem.M = M;
        problem.name = 'bio8';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        optimum = W1_Optimum(500,M);
        problem.fnc = @(x)bio8(x,problem.M,optimum);
        problem.F_best = optimum;
    case 9 
        dim = d;
        problem.M = M;
        problem.name = 'bio9';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        optimum = W2_Optimum(500,M);
        problem.fnc = @(x)bio9(x,problem.M,optimum);
        problem.F_best = optimum;
    case 10 
        dim = d;
        problem.M = M;
        problem.name = 'bio10';
        problem.dims = dim;
        problem.lb = 0.1 * all_kcat_values;
        problem.ub = 1.0 * all_kcat_values;
        optimum = W3_Optimum(500,M);
        problem.fnc = @(x)bio10(x,problem.M,optimum);
        problem.F_best = optimum;
        
%     case 11
%         dim = d;
%         problem.M = M;
%         problem.name = 'WFG4';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = W4_Optimum(500,M);
%         problem.fnc = @(x)WFG4(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 12 
%         dim = d;
%         problem.M = M;
%         problem.name = 'WFG5';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = W5_Optimum(500,M);
%         problem.fnc = @(x)WFG5(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 13 
%         dim = d;
%         problem.M = M;
%         problem.name = 'ZDT1';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = Z1_Optimum(500,M);
%         problem.fnc = @(x)ZDT1(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 14 
%         dim = d;
%         problem.M = M;
%         problem.name = 'ZDT2';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = Z2_Optimum(500,M);
%         problem.fnc = @(x)ZDT2(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 15 
%         dim = d;
%         problem.M = M;
%         problem.name = 'ZDT3';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = Z3_Optimum(500,M);
%         problem.fnc = @(x)ZDT3(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 16 
%         dim = d;
%         problem.M = M;
%         problem.name = 'ZDT4';
%         problem.dims = dim;
%         problem.lb = [0, -5*ones(1, dim-1)];
%         problem.ub = [1,  5*ones(1, dim-1)];
%         optimum = Z4_Optimum(500,M);
%         problem.fnc = @(x)ZDT4(x,problem.M,optimum);
%         problem.F_best = optimum;
%     case 17 
%         dim = d;
%         problem.M = M;
%         problem.name = 'ZDT6';
%         problem.dims = dim;
%         problem.lb = 0*ones(1,dim);
%         problem.ub = 1*ones(1,dim);
%         optimum = Z6_Optimum(500,M);
%         problem.fnc = @(x)ZDT6(x,problem.M,optimum);
%         problem.F_best = optimum;
   
end
