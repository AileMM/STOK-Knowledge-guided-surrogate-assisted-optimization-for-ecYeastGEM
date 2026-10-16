function [best_candidate,best_candidate_obj,paras] =transfer_opt(database,name_saea,paras,func_surrogate,M,transfer_source,knowledge_base)
% 设置不同优化器，框架式改为不同的算法
        most_similar = [];
switch name_saea
    % NSGAII
    case 'NSGAII'
        %% 构建surrogate的时候需要计算相似度开始进行迁移
        configs.database = database; %database的实际值
        configs.popsize = 50;
        configs.query_max = 50;
        configs.ini_strategy = 'elite+random';
        configs.acquisition = 'lcb';
        configs.search_engine = 'ea';
        configs.selector = 'roulette_wheel'; 
      
       %% 种群更新--产生子代用代理预测 模型预测值本就是LCB值
        [population,objs_acq] = transfer_pre(func_surrogate,configs,M,transfer_source);
        solution_promising = population;
        objs_promising = objs_acq;
  
   
       %% 初始化一个map划分源任务类型
        grouped_KB = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for i = 1:length(knowledge_base)
            mm = size(knowledge_base(i).objs, 2);  % 当前任务的目标维度
            key = num2str(mm);                     % 用目标维度作为 key，例如 '2', '3', ...
            if isKey(grouped_KB, key)
                tmp = grouped_KB(key);     % 先取出来（是个 cell）
                tmp{end+1} = knowledge_base(i);  % 更新 cell
                grouped_KB(key) = tmp;     % 再放回去
            else
                grouped_KB(key) = {knowledge_base(i)};
            end
        end
        M_query = M;
        key = num2str(M_query);
        source_tasks = grouped_KB(key);  

        [sorted_idx, ~] = compare_PFs_by_cosine(objs_promising, source_tasks);
        most_similar(end+1) = sorted_idx;   % 记录这次的索引
        disp(most_similar);
        
        source_PF = source_tasks{sorted_idx}.objs;
        % 判断哪些预测个体不被 source_PF 中任何一个个体支配
        n = size(objs_promising, 1);
        is_candidate = true(n, 1);
        for i = 1:n
            for j = 1:size(source_PF, 1)
                if dominates(source_PF(j, :), objs_promising(i, :))
                    is_candidate(i) = false;
                    break;
                end
            end
        end
        % 提取未被支配的个体作为候选
        selected_idx = find(is_candidate);
        best_candidate = solution_promising(selected_idx, :);
        best_candidate_obj = objs_promising(selected_idx, :);
%         
   
end
end
function flag = dominates(a, b)
    flag = all(a <= b) && any(a < b);  % a 支配 b
end