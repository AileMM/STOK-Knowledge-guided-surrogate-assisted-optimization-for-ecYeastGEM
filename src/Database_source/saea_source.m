function [solutions,objs, surrogate, selected_variables, selected_C, selected_D]= saea_source(problem,FEsMax,M,id_source)

% saea parameters
popsize = 50; % the population size
% surrogate = 'gpr'; % the surrogate model
% surrogate = 'prs'; % the surrogate model
% surrogate = 'rbf'; % the surrogate model
% query_max = 50; % =1: preselection, >1: iteration
% selector = 'roulette_wheel'; % the selestion mechanism
acquisition = 'lcb'; % the lower confidence bound for acquisition
ini_strategy = 'elite+random'; % the strategy for initializing the population at each round of surrogate-assisted search
% search_engine = 'ea';
num_initial_solutions = 50; % the number of initial solutions

% initialization
fun = problem.fnc; % the objective function
lb = problem.lb; % the lower bound
ub = problem.ub; % the upper bound
solutions_initial = lhsdesign_modified(num_initial_solutions,lb,ub); %% 改成生物模型初始化，按照变量缩放修改

% 这里需要传递一个M目标个数初始化目标值
objs_initial = zeros(num_initial_solutions,M);
for i = 1:num_initial_solutions
    objs_initial(i,:) = fun(solutions_initial(i,:));  %计算初始目标值
end
FEsUsed = num_initial_solutions;
database = [solutions_initial, objs_initial];
flag = 0;
X = database(:,1:end-M);  %需要传递一个M获得X和Y
y = database(:,end-M+1:end);
while FEsUsed < FEsMax
    %% 训练不同代理模型
    while flag == 0  
        try
             % 预分配存储不同的代理模型
            func_surrogate = struct();
            % 训练GPR
            func_surrogate.gpr = cell(1, M);
            for i = 1:M
                func_surrogate.gpr{i} = surrogate_model(X, y(:, i), 'gpr', M);
            end
            % 训练kri
            func_surrogate.kri = cell(1, M);
            for i = 1:M
                func_surrogate.kri{i} = surrogate_model(X, y(:, i), 'kriging', M);
            end
            % 训练RBF
            func_surrogate.rbf = surrogate_model(X, y, 'rbf', M);
            % 设置 flag，表示训练成功，退出 while 循环
            flag = 1;
        catch
            continue;
        end
    end  
    % identify a candidate solution to be evaluated by applying the evolutionary solver to the surrogate model
    % 现在 objs_acq，objs_acq_offspring都是NxMx3的数据，包含有三个代理模型的预测结果
    [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,func_surrogate,acquisition,M);
    rmse = zeros(size(objs_acq, 3),1);
    for i =1:size(objs_acq, 3)
        temp = objs_acq(:, :, i);  % 先存储索引结果
        rmse(i) = sqrt(mean((y(:) - temp(:)).^2));  % 再进行计算
    end
    % =======计算RMSE决定使用哪个代理模型对子代预测========
    [~,idx] = min(rmse);
    fields = fieldnames(func_surrogate);  
    values = struct2cell(func_surrogate);
    nth_field_name = fields{idx}; 
    nth_field_value = values{idx};
    new_surrogate.(nth_field_name) = nth_field_value;
    % ======================================================
    population_offspring = Generation(population); 
    objs_acq_offspring = sur_Evaluation(new_surrogate,population_offspring,acquisition,M);  %模型预测
    Dec = [population; population_offspring];
    Obj = [y; objs_acq_offspring];
    [rank,FrontNo] = NSGASort(Obj);
    population = Dec(rank(1:num_initial_solutions),:);
    objs_acq = Obj(rank(1:num_initial_solutions),:);
    
    
%     % evaluate the candidate solution and update the database
%     epsilon = 5e-3;
%     no_trials  = 50;
%     scales = linspace(0.1,1,no_trials);
%     c = 1;
%     while min(max(abs(repmat(obj_candidate, size(database,1), 1) - database(:, end-M+1:end)), [], 2)) < epsilon
%         % 计算扰动=======这里后期实现目的后再进行扰动吧，决策空间和目标空间调整一下
%         perturbation = scales((mod(c,no_trials)==0)*no_trials+(mod(c,no_trials)~=0)*mod(c,no_trials)) ...
%                        * (ub-lb) .* (rand(1,dim) - 0.5);
%         % 施加扰动
%         solution_candidate = solution_candidate + perturbation;
%         % 确保新解仍在边界范围内
%         solution_candidate = max(min(solution_candidate, ub), lb);
%         % 计算新解在目标空间的值（假设 `evaluate` 是目标函数计算方法）
%         obj_candidate = evaluate(solution_candidate);  
%         % 计数器更新
%         c = c + 1;
%     end
%     obj_candidate = fun(solution_candidate);
%     database = [database;solution_candidate obj_candidate];
    [~,FrontNo] = NSGASort(objs_acq);
    idx = FrontNo==1;
    solution_candidate = population(idx,:);
    obj_candidate = fun(solution_candidate); 
    new_data = [solution_candidate, obj_candidate];
    new_database = [database; new_data];
    [solutions,objs] = Entropy(new_database, M);
%     solutions = (database(:,1:end-M) - repmat(lb, size(database,1), 1)) ./ (repmat(ub, size(database,1), 1) - repmat(lb, size(database,1), 1));
%     objs = database(:, end-M+1:end);
    surrogate = new_surrogate;
    FEsUsed = FEsUsed+ size(solution_candidate,1);
end


% 数值梯度（Numerical Gradient）前向差分法（Forward Finite Difference Method）
X_nd = solutions; % 决策变量
Y_nd = objs; % 目标函数值
% 计算梯度（每个决策变量对目标的影响）
grad = zeros(size(X_nd,1), size(X_nd,2), size(Y_nd,2));  % N × D × M 维度
delta = 1e-5;
for i = 1:size(X_nd,2)  % 遍历每个决策变量
    X_perturbed = X_nd;
    X_perturbed(:,i) = X_nd(:,i) + delta;  % 扰动第i维
%     deltaY = fun(X_perturbed) - fun(X_nd);  % N × M
    grad(:,i,:) = (fun(X_perturbed) - fun(X_nd)) / delta;  % 计算每个目标的梯度
end
grad_mean = squeeze(mean(abs(grad), [1, 3]));  % D×1 所有目标的平均梯度

% Sobol敏感度分析
N = size(X_nd, 1);  % 采样点数
D = size(X_nd, 2);  % 变量维度
score = zeros(1, D); 
for i = 1:D
    X_permuted = X_nd;  
    X_permuted(:, i) = X_permuted(randperm(N), i);  % 随机打乱第 i 维
    Y_permuted = fun(X_permuted);  
    a=var(Y_permuted);
    b=var(Y_nd);
    disp(a);
    disp(b);
    delta_var  = abs(var(Y_permuted) - var(Y_nd)) ./ var(Y_nd);
    score(i) = mean(delta_var);
end

% 整体输出变量分析 根据阈值判断
score = score / max(score);
grad_score  = grad_mean / max(grad_mean); 
final_score = 0.4 * grad_score  + 0.6 * score;  
selected_variables = find(final_score > 0.2);  

%% 到收敛和多样变量分类
% 梯度方向一致性（低 std 代表一致）
grad_sign = sign(grad);  % 计算梯度的符号
dir_std = squeeze(std(grad_sign, 0, 3));   % N × D，计算每个变量的梯度方向的标准差
dir_consistency = 1 - mean(dir_std, 1);    % 越大越一致，取 1 - std 得到一致性度量
selected_converge = find(dir_consistency > 0.5); 

% 多样性分析
grad_std = squeeze(std(grad, 0, 1));         % D × M
diversity_score = mean(abs(grad_std), 2);    % D × 1
diversity_score = diversity_score / max(diversity_score);
selected_diversity = (find(diversity_score > 0.5))';

% 输出
selected_C = intersect(selected_variables, selected_converge);
selected_D = intersect(selected_variables, selected_diversity);
% selected_vars = unique([selected_variables, selected_C, selected_D]);








