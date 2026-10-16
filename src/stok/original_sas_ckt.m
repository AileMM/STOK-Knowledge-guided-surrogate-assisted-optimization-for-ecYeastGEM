function [Out_dec,Out_obj, metric, paras]= sas_ckt(fun,lb,ub,FEsMax,name_saea,knowledge_base,paras,M,optimum)

    % =================Phase 1: Initialization===================
    if isempty(name_saea)
        name_saea = 'NSGAII';
    end
    num_sources = length(knowledge_base); 
    interval_transfer = paras.interval;
    num_initial_solutions = size(paras.pop,1); % 初始化种群个数
    solutions_initial = paras.pop;   % 初始种群
    objs_initial = zeros(num_initial_solutions,M);
    %% 如果后面变量需要迁移的话，这里就不需要计算目标函数值
    for i = 1:num_initial_solutions
        objs_initial(i,:) = fun(solutions_initial(i,:));  %fun是传递的目标任务
    end
    FEsUsed = num_initial_solutions;
    database = [solutions_initial, objs_initial]; % 有目标值的种群
%     dim = length(lb);
    transfer_states = [];
    if strcmp(name_saea,'AutoSAEA')
        load('optimizers\auto-saea\autosaea_paras.mat');
        options_saea.func_exp = fun;
    else
        options_saea = [];
    end

%     % build the surrogate models of k source tasks
%     %% 为每个源任务构建高斯过程代理模型
%     surrogates_source = struct; 
%     n_count = 0;
%     while n_count < num_sources
%         i = n_count + 1;
%         try 
%             surrogates_source(i).func = knowledge_base(i).surr;
%             n_count = n_count+1;
%         catch
%             continue;
%         end
%     end
%     %% execute the task adaptation if needed 这一步为了解决任务异构
%     %如果 ada_flag 为 1，则通过 solution_adaptation 方法计算任务之间的迁移矢量（ada_vectors）
%     %% 迁移矢量用于将源任务知识适配到目标任务。
%     if paras.ada_flag == 0
%         ada_vectors = zeros(num_sources,dim);
%     else
%         solutions_target = (solutions_initial-repmat(lb,num_initial_solutions,1))./...
%                 (repmat(ub,num_initial_solutions,1)-repmat(lb,num_initial_solutions,1));
%         objs_target = objs_initial;
%         ada_vectors = solution_adaptation(surrogates_source,solutions_target,objs_target);
%     end

    Out_dec = database(:,1:end-M); 
    Out_obj = database(:,end-M+1:end);
    
    %% 设置一个存档保存非支配解作为最后的输出
    [~,FrontNo] = NSGASort(Out_obj);
    index = FrontNo==1;
    archive_Dec = Out_dec(index,:);
    archive_obj = Out_obj(index,:);
    
    %% 初始化存档参数设置
    best_IGD = inf;
    phase_mode = 'self'; 
    igd_last = inf;  % 初始化前一代 IGD
    isFirstWarmStart = true;  % 标志是否第一次进入
    
    window_size = 10; 
    progress_threshold = 1e-3;
    %% 连续几代 IGD 未改善时触发迁移
    stall_threshold = 50;  
    stall_counter = 0;  

    transfer_fail_counter = 0;         % 迁移连续失败的次数
    transfer_fail_threshold = 10;      % 达到这个就切回 self opt

    igd_no_change_counter = 0;         % IGD连续无提升的计数器
    igd_no_change_threshold = 200;

    epsilon_igd = 1e-3;               % 认为有“明显”IGD提升的标准
    
    IGD_record = [];
    FE_record = [];
    IGD_history=[];
    max_transfer_generations = 5;       
    archive_best_obj = archive_obj;
    archive_best_dec = archive_Dec;
    archive_limit = 500; 

    % ==========================Phase 2: Optimization=====================
    while FEsUsed < FEsMax
%         % 测试集
%         pop = rand(20, 10);
%         pop_obj = fun(pop);
%         
%         for i = 1:M
%             func_surrogate.gpr{i} = surrogate_model(Out_dec, Out_obj(:, i), 'gpr', M);
%         end
%         % 计算预测误差
%         pred_error1 = zeros(1, M);
%         for i = 1:M
%             Y_pred = func_surrogate.gpr{i}.func(pop);  % GPR模型的预测函数
%             pred_error1(i) = mean((Y_pred - pop_obj(:, i)).^2); 
%         end
%         mean_pred_error1 = mean(pred_error1); 
%         disp(mean_pred_error1);

%         [newdatabase,FEs_inc,success,transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
%         if success
%             FEsUsed = FEsUsed + FEs_inc;
%         end
        
       % Step 1: 训练/更新代理模型 
         if isFirstWarmStart
            % 第一次冷启动，只选 transfer source，不迁移个体为什么？？
            [newdatabase, FEs_inc, ~, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
            database = newdatabase;
            FEsUsed = FEsUsed + FEs_inc;
            isFirstWarmStart = false;
        else
            target_pop = database(:,1:end-M);
            objs = database(:,end-M+1:end);
            func_surrogate.gpr = cell(1, M);
            for i = 1:M
                func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
            end
            % Step 2: 评估代理模型误差
            current_rmse = evaluate_model_error(func_surrogate, database,M);
            % Step 3: 判断是否允许迁移
            threshold_rmse = 0.05;  % 可以根据实验调整
            if current_rmse > threshold_rmse
                % 代理模型还不够准，继续warm start
                [newdatabase, FEs_inc, success, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
                database = newdatabase;
                FEsUsed = FEsUsed + FEs_inc;
                if success
                    update_model_flag = true;
                end
            else
                % 代理模型已经够准了，禁止warm start
                success = false;
            end
         end
         
        target_pop = newdatabase(:,1:end-M);
        objs = newdatabase(:,end-M+1:end);
        func_surrogate.gpr = cell(1, M);
        for i = 1:M
            func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
        end
        
%         % 计算预测误差
%         pred_error2 = zeros(1, M);
%         for i = 1:M
%             Y_pred = func_surrogate.gpr{i}.func(pop);  % GPR模型的预测函数
%             pred_error2(i) = mean((Y_pred - pop_obj(:, i)).^2);  % MSE
%         end
%         mean_pred_error2 = mean(pred_error2);  % 所有目标的平均误差
%         disp(mean_pred_error2);
        %% 自迁移
        [solution_in,improvement_in,options_saea] = self_opt(database,name_saea,options_saea,func_surrogate,M,archive_best_obj);
        %% 动态迁移强度设置
        if strcmp(phase_mode, 'transfer')
            [solution_ex,improvement_ex] = transfer_opt(database,name_saea,options_saea,func_surrogate,M,transfer_source,knowledge_base);
            relation = dominance_relation(improvement_in,improvement_ex);
            if relation == 1
                solution_candidate = solution_in;
                solution_candidate_obj = improvement_in;
            elseif relation == -1
                solution_candidate = solution_ex;
                solution_candidate_obj = improvement_ex;
            else
                solution_candidate = [solution_in;solution_ex];
                solution_candidate_obj = [improvement_in; improvement_ex];
%                 solution_candidate = solution_ex;
%                 solution_candidate_obj = improvement_ex;
            end
        else
            solution_candidate = solution_in;
            solution_candidate_obj = improvement_in;
            disp('Current candidate solution:');
            disp(solution_candidate);
        end
%         if mod(FEsUsed,interval_transfer)==0 
%         if stall_counter >= stall_threshold
%             [solution_ex,improvement_ex] = transfer_opt(newdatabase,name_saea,options_saea,func_surrogate,M,transfer_source,knowledge_base);
%             relation = dominance_relation(improvement_in,improvement_ex);
%             if relation == 1
%                 solution_candidate = solution_in;
%                 solution_candidate_obj = improvement_in;
%             elseif relation == -1
%                 solution_candidate = solution_ex;
%                 solution_candidate_obj = improvement_ex;
%             else
%                 solution_candidate = [solution_in;solution_ex];
%                 solution_candidate_obj = [improvement_in; improvement_ex];
%             end
%         else
%             solution_candidate = solution_in;
%             solution_candidate_obj = improvement_in;
%             disp('Current candidate solution:');
%             disp(solution_candidate);
%         end

        if strcmp(phase_mode, 'self')
            % 这是一个扰动
            solution_candidate = diversify_candidate(solution_candidate, archive_best_dec, lb, ub);
        end
% 
%         % 2. 判断哪些个体未被 archive_best_obj 支配
%         mask = dominance_relation(solution_candidate_obj, archive_best_obj);  % 上面定义过的函数
%         % 3. 只保留未被支配的解，进行真实函数评估
%         solution_candidate_filtered = solution_candidate(mask, :);
        
        solution_candidate_filtered = solution_candidate;
        
       %% 设置迁移数目
        K = 8;
        num = min(K, size(solution_candidate_filtered,1));
        selected_idx = randperm(size(solution_candidate_filtered,1), num);
        solution_candidate_filtered = solution_candidate_filtered(selected_idx, :);
        if ~isempty(solution_candidate_filtered)
            obj_candidate = fun(solution_candidate_filtered);
            Out_temp_dec = [Out_dec; solution_candidate_filtered]; 
            Out_temp_obj = [Out_obj; obj_candidate];
            TEMP = [Out_temp_dec, Out_temp_obj];
            FEsUsed = FEsUsed + size(solution_candidate_filtered,1);
        else
            Out_temp_dec = Out_dec;
            Out_temp_obj = Out_obj;
            TEMP = [Out_temp_dec, Out_temp_obj];
        end
        
        % =======种群更新======
        [rank,~,~] = NSGASort(Out_temp_obj);
        database = TEMP(rank(1:num_initial_solutions),:);
        
        % ======更新archive=archive+population======
        combined_dec = [archive_best_dec; Out_temp_dec];
        combined_obj = [archive_best_obj; Out_temp_obj];

        % 去除支配解，保留真正的非支配个体
        [~, FrontNo_combined] = NSGASort(combined_obj);
        idx_best = FrontNo_combined == 1;
        % 提取合并后的非支配解
        archive_best_dec = combined_dec(idx_best, :);
        archive_best_obj = combined_obj(idx_best, :);
     
%         if size(archive_best_obj,1) > archive_limit
%             crowd_dis = CrowdingDistance(archive_best_obj);
%             [~, sort_idx] = sort(crowd_dis, 'descend');
%             selected = sort_idx(1:archive_limit);
%             archive_best_dec = archive_best_dec(selected, :);
%             archive_best_obj = archive_best_obj(selected, :);
%         else
%             archive_best_dec = combined_dec;
%             archive_best_obj = combined_obj;
%         end

       % ===========计算指标============
        IGD = calIGD(archive_best_obj, optimum);
        if IGD < best_IGD
            best_IGD = IGD;
        end
        IGD_history = [IGD_history, IGD];
%         % ==== 阶段控制逻辑 ====
%         if strcmp(phase_mode, 'self')
%             if abs(IGD - igd_last) < 1e-6
%                 igd_no_change_counter = igd_no_change_counter + 1;
%             else
%                 igd_no_change_counter = 0;
%             end
%             igd_last = IGD;
%             stall_counter = stall_counter + 1;
%             
% %             if igd_no_change_counter >= igd_no_change_threshold || stall_counter >= stall_threshold
%             if  stall_counter >= stall_threshold
%                 fprintf('>>> Switch to transfer phase...\n');
%                 phase_mode = 'transfer';
%                 transfer_fail_counter = 0;
%                 last_transfer_IGD = IGD;
%             end
%         elseif strcmp(phase_mode, 'transfer')
% %             if IGD < last_transfer_IGD - epsilon_igd
%             if IGD < last_transfer_IGD
%                 fprintf('>>> Transfer success! Stay in transfer phase.\n');
%                 last_transfer_IGD = IGD;
%                 transfer_fail_counter = 0;
%                 phase_mode = 'self';   % 这里一旦成功了，就直接转为自优化模式
%             else
%                 transfer_fail_counter = transfer_fail_counter + 1;
%                 fprintf('>>> Transfer fail (%d/%d)\n', transfer_fail_counter, transfer_fail_threshold);
%                 if transfer_fail_counter >= transfer_fail_threshold
%                     fprintf('>>> Exit transfer. Back to self opt.\n');
%                     phase_mode = 'self';
%                     stall_counter = 0;
%                     igd_no_change_counter = 0;
%                 end
%             end
%         end
%         if mod(current_generation, window_size) == 0
%             recent_progress = IGD_history(end-window_size+1) - IGD_history(end);
%             progress_rate = recent_progress / IGD_history(end-window_size+1);
% 
%             if phase_mode == "self"
%                 if progress_rate < progress_threshold
%                     phase_mode = "transfer";
%                     transfer_counter = 0;
%                     fprintf(">>> Enter transfer phase due to low progress\n");
%                 end
%             elseif phase_mode == "transfer"
%                 transfer_counter = transfer_counter + 1;
%                 if transfer_counter >= max_transfer_generations
%                     phase_mode = "self";
%                     fprintf(">>> Return to self-optimization after transfer window\n");
%                 end
%             end
%         end
        if length(IGD_history) >= window_size && mod(FEsUsed, window_size) == 0
            idx_start = length(IGD_history) - window_size + 1;
            recent_progress = IGD_history(idx_start) - IGD_history(end);
            progress_rate = recent_progress / (IGD_history(idx_start) + eps);  % 加eps防除0

            if phase_mode == "self"
                if progress_rate < progress_threshold
                    phase_mode = "transfer";
                    transfer_counter = 0;
                    fprintf(">>> Enter transfer phase due to low progress\n");
                end
            elseif phase_mode == "transfer"
                transfer_counter = transfer_counter + 1;
                if transfer_counter >= max_transfer_generations
                    phase_mode = "self";
                    fprintf(">>> Return to self-optimization after transfer window\n");
                end
            end
        end


        IGD_record = [IGD_record; best_IGD];
        FE_record = [FE_record; FEsUsed];
        if paras.ada_flag == 0
            fprintf([name_saea,'-CKT-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,best_IGD);
        else
            fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,IGD);
        end
    end
    % =====输出archive=======
    Out_dec = archive_best_dec;
    Out_obj = archive_best_obj;
    metric = [FE_record,IGD_record];
end

function dis = CrowdingDistance(obj)
    [N, M] = size(obj);
    dis = zeros(N,1);
    Fmax = max(obj,[],1);
    Fmin = min(obj,[],1);
    for i = 1:M
        [~, rank] = sort(obj(:,i));
        dis(rank(1)) = inf;
        dis(rank(end)) = inf;
        for j = 2:N-1
            if Fmax(i) - Fmin(i) == 0
                d = 0;
            else
                d = (obj(rank(j+1),i) - obj(rank(j-1),i)) / (Fmax(i) - Fmin(i));
            end
            dis(rank(j)) = dis(rank(j)) + d;
        end
    end
end


function result = dominance_relation(A, B)
% A: n × M matrix (多个候选解)
% B: m × M matrix (源任务的 PF)
% result: n × 1 logical 向量，true 表示 A 中对应解未被 B 中任何一个支配

    n = size(A, 1);
    m = size(B, 1);
    result = true(n, 1);  % 初始化所有解为未被支配

    for i = 1:n
        for j = 1:m
            if all(B(j, :) <= A(i, :)) && any(B(j, :) < A(i, :))
                result(i) = false;  % 一旦被支配，标记为 false
                break;  % 不再检查其他 j
            end
        end
    end
end
