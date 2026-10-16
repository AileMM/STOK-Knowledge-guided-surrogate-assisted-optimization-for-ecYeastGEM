function [Out_dec,Out_obj, metric, all_avg_norm, reachedFE,runtime,fidelity]= sas_ckt(fun,lb,ub,FEsMax,name_saea,knowledge_base,paras,M,optimum,name)

    % Paired, pre-evaluation surrogate fidelity records (optional output #7).
    % Columns are objectives in the same units/order returned by fun.
    fidelity = struct();
    fidelity.candidate = zeros(0, size(paras.pop,2));
    fidelity.predicted = zeros(0,M);
    fidelity.selection_obj = zeros(0,M); % self_opt/transfer_opt output, unmodified
    fidelity.exact = zeros(0,M);
    fidelity.iteration = zeros(0,1);
    fidelity.phase = cell(0,1);
    fidelity.selected_idx = zeros(0,1);
    fidelity.batch_row = zeros(0,1);
    fidelity.FE_before = zeros(0,1);
    fidelity.FE = zeros(0,1); % logical candidate FE index, not solver-call count
    fidelity.FE_after = zeros(0,1);
    fidelity.training_size = zeros(0,1);
    fidelity.finite_pair = false(0,M);
    % ID mapping for the candidate knowledge base passed by main.
    fidelity.target_id = [];
    fidelity.source_ids = [];
    if isfield(paras, 'target_id')
        fidelity.target_id = paras.target_id;
    end
    if isfield(paras, 'source_ids')
        fidelity.source_ids = paras.source_ids;
    end
    fidelity.variant = 'Full';
    fidelity.surrogate_transfer_enabled = true;
    fidelity.transfer_opt_enabled = true;
    fidelity.task = name;
    fidelity.optimizer = name_saea;
    fidelity.initial_rng = rng;
    fidelity.FE_max = FEsMax;
    fidelity.completed = false;
    fidelity_file = '';
    if isfield(paras, 'fidelity_file') && ~isempty(paras.fidelity_file)
        fidelity_file = paras.fidelity_file; % use a distinct path per run
    end
    iteration = 0;

    % =================Phase 1: Initialization===================
    if isempty(name_saea)
        name_saea = 'NSGAII';
    end
%     num_sources = length(knowledge_base); 
%     interval_transfer = paras.interval;
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

    Out_dec = database(:,1:end-M); 
    Out_obj = database(:,end-M+1:end);
    
    %% 设置一个存档保存非支配解作为最后的输出
    [~,FrontNo] = NSGASort(Out_obj);
    index = FrontNo==1;
    archive_Dec = Out_dec(index,:);
    archive_obj = Out_obj(index,:);
    
    IGD = calIGD(archive_obj, optimum);
    fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,IGD);
    
    %% 初始化存档参数设置
    best_IGD = inf;
    phase_mode = 'self'; 
    igd_last = inf;  % 初始化前一代 IGD
    isFirstWarmStart = true;  % 标志是否第一次进入
    
    norm_flag = false;
    norm_obj_history = [];
    norm_obj_record = [];
    all_avg_norm = [];
    
    %% 迁移转换初始化    
    igd_history = [];
    window_size = 10;
    improvement_threshold = 1e-4;
    max_transfer_attempts = 5;
    transfer_attempts = 0;
    cooldown_FE = 50;        % 冷却阶段，单位是函数评估次数（可调）
    last_switch_FE = 0;       % 上次阶段切换的 FEs
    
    IGD_record = [];
    FE_record = []; 
    archive_best_obj = archive_obj;
    archive_best_dec = archive_Dec;
    archive_limit = 500; 
    
    % 加速比求解
    reachedFE = [];
    threshold_map = containers.Map( ...
    {'bio1', 'bio2', 'bio3', 'bio4','bio5','bio6','bio7','bio8','bio9','bio10'}, ...
    {[0.3, 3], [0.3, 4], [0.32, 4], [0.29, 10], [0.36, 7], [0.34, 5], [0.35, 2], [0.3, 25], [0.3, 10], [0.3, 11]} ...
);
    
%     %% 记录模型的性能提升初始化集合
%     rmseA_record = [];
%     rmseB_record = [];
%     maeA_record = [];
%     maeB_record = [];
%     FE_Record = [];
%     R2A_record = []; 
%     R2B_record = [];

    % =========================Phase 2: Optimization=====================
%     %% 测试集生成验证模型性能
%     Ytest = zeros(num_initial_solutions,M);
%     num_samples = 100;      % 测试集样本数量
%     num_vars = 20;          % 决策变量维度
%     lower_bounds = zeros(1, num_vars); 
%     upper_bounds = ones(1, num_vars);   
%     % Latin Hypercube Sampling
%     lhs_unit = lhsdesign(num_samples, num_vars);  % 在 [0,1]^n 中采样
%     X_test = repmat(lower_bounds, num_samples, 1) + ...
%      lhs_unit .* repmat((upper_bounds - lower_bounds), num_samples, 1);
%      for i = 1:num_samples
%         Ytest(i,:) = fun(X_test(i,:)); 
%      end
    tic;
    while FEsUsed < FEsMax
        iteration = iteration + 1;
%          % ------未迁移之前的数据集------
%         target_pop = database(:,1:end-M);
%         objs = database(:,end-M+1:end);
%         func_surrogateA.gpr = cell(1, M);
%         for i = 1:M
%             func_surrogateA.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
%         end
%         for i = 1:M
%             [objs_pre(:, i),objs_std(:, i)] = func_surrogateA.gpr{i}.func(X_test);
%         end
%         Y_predA = objs_pre;
       % Step 1: 训练/更新代理模型 
         if isFirstWarmStart
            [newdatabase, FEs_inc, ~, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
            % [newdatabase, FEs_inc, ~, transfer_source] = warm_start_variable_transfer_LAE(database, knowledge_base,fun, M);
            
            
            database = newdatabase;
            FEsUsed = FEsUsed + FEs_inc;
            isFirstWarmStart = false;
         end    
        target_pop = database(:,1:end-M);
        objs = database(:,end-M+1:end);
        func_surrogate.gpr = cell(1, M);
        for i = 1:M
            func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
        end
        fidelity_training_size = size(target_pop,1);
        if ~isFirstWarmStart
            % Step 2: 评估代理模型误差
            current_rmse = evaluate_model_error(func_surrogate, database,M);
            disp(current_rmse);
            %% Step 3: 判断是否允许迁移 代理模型的精度这里需要消融
%             threshold_rmse = 0.05;  
            threshold_rmse = 1e-4*5; 
            if current_rmse > threshold_rmse
                % 代理模型还不够准，继续warm start
                [newdatabase, FEs_inc, success, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
                database = newdatabase;
                FEsUsed = FEsUsed + FEs_inc;
                if success
                    fidelity_training_size = size(database,1);
                    update_model_flag = true;
                    % warm start 后更新 surrogate 模型
                    target_pop = database(:,1:end-M);
                    objs = database(:,end-M+1:end);
                    for i = 1:M
                        func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
                    end
                end
            end
        end
%         %% =====计算 RMSE R2 MAE=====
%         target_pop = database(:,1:end-M);
%         objs = database(:,end-M+1:end);
%         func_surrogateB.gpr = cell(1, M);
%         for i = 1:M
%             func_surrogateB.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
%         end
%         for i = 1:M
%             [objs_pre(:, i),objs_std(:, i)] = func_surrogateB.gpr{i}.func(X_test);
%         end
%         Y_predB = objs_pre;
%         rmseA = sqrt(mean((Y_predA - Ytest).^2));
%         rmseB = sqrt(mean((Y_predB - Ytest).^2));
%         maeA = mean(abs(Y_predA - Ytest));
%         maeB = mean(abs(Y_predB - Ytest));
%         % R2 for Model A
%         SS_res_A = sum((Ytest - Y_predA).^2);                  % Residual sum of squares
%         SS_tot_A = sum((Ytest - mean(Ytest)).^2);              % Total sum of squares
%         R2_A = 1 - SS_res_A ./ SS_tot_A;                       % R^2 per objective
%         % R2 for Model B
%         SS_res_B = sum((Ytest - Y_predB).^2);
%         SS_tot_B = sum((Ytest - mean(Ytest)).^2);
%         R2_B = 1 - SS_res_B ./ SS_tot_B;
%         % 记录
%         rmseA_record(end+1,:) = rmseA;
%         rmseB_record(end+1,:) = rmseB;
%         maeA_record(end+1,:) = maeA;
%         maeB_record(end+1,:) = maeB;
%         R2A_record(end+1, :) = R2_A;
%         R2B_record(end+1, :) = R2_B;
%         FE_Record(end+1) = FEsUsed;
%         % 打印
%         fprintf('FEs = %d | RMSE A: %.4f, B: %.4f | MAE A: %.4f, B: %.4f\n', ...
%             FEsUsed, mean(rmseA), mean(rmseB), mean(maeA), mean(maeB));

        %% 自迁移
        [solution_in,improvement_in,options_saea] = self_opt(database,name_saea,options_saea,func_surrogate,M,archive_best_obj);
%         solution_candidate = solution_in;
%         solution_candidate_obj = improvement_in;
        %% 动态迁移强度设置
        if strcmp(phase_mode, 'transfer')
            [solution_ex,improvement_ex] = transfer_opt(database,name_saea,options_saea,func_surrogate,M,transfer_source,knowledge_base);
            idx_B = filter_nondominated(improvement_ex, improvement_in);
            solution_candidate = solution_ex(idx_B, :);
            solution_candidate_obj = improvement_ex(idx_B, :);
        else
            solution_candidate = solution_in;
            solution_candidate_obj = improvement_in;
        end
        solution_candidate_filtered = solution_candidate;

       %% 填充个体数
        K = 11;
        num = min(K, size(solution_candidate_filtered,1));
        selected_idx = randperm(size(solution_candidate_filtered,1), num);
        solution_candidate_filtered = solution_candidate_filtered(selected_idx, :);
        if ~isempty(solution_candidate_filtered)
            n_candidate = size(solution_candidate_filtered,1);
            assert(isequal(size(solution_candidate_obj), ...
                [size(solution_candidate,1), M]), ...
                'STOK:SelectionShape', 'Candidate selection objectives must be N-by-M.');
            selection_obj_filtered = solution_candidate_obj(selected_idx,:);
            % Evaluate the frozen current model BEFORE exact evaluation/retraining.
            % Do not assume improvement_in/ex are raw surrogate objective means.
            predicted_candidate = zeros(n_candidate,M);
            for objective_idx = 1:M
                prediction = func_surrogate.gpr{objective_idx}.func(solution_candidate_filtered);
                assert(isvector(prediction) && numel(prediction) == n_candidate, ...
                    'STOK:PredictionShape', 'Surrogate must return one value per candidate.');
                predicted_candidate(:,objective_idx) = prediction(:);
            end
            FE_before_candidate = FEsUsed;
            obj_candidate = fun(solution_candidate_filtered);
            assert(isequal(size(obj_candidate), [n_candidate,M]), ...
                'STOK:ExactShape', 'fun must return an N-by-M matrix in input row order.');
            % Append before sorting/truncation; retain even nonfinite pairs for audit.
            rows = size(fidelity.exact,1) + (1:n_candidate);
            fidelity.candidate(rows,:) = solution_candidate_filtered;
            fidelity.predicted(rows,:) = predicted_candidate;
            fidelity.selection_obj(rows,:) = selection_obj_filtered;
            fidelity.exact(rows,:) = obj_candidate;
            fidelity.iteration(rows,1) = iteration;
            fidelity.phase(rows,1) = repmat({phase_mode},n_candidate,1);
            fidelity.selected_idx(rows,1) = selected_idx(:);
            fidelity.batch_row(rows,1) = (1:n_candidate)';
            fidelity.FE_before(rows,1) = FE_before_candidate;
            fidelity.FE(rows,1) = FE_before_candidate + (1:n_candidate)';
            fidelity.FE_after(rows,1) = FE_before_candidate + n_candidate;
            fidelity.training_size(rows,1) = fidelity_training_size;
            fidelity.finite_pair(rows,:) = isfinite(predicted_candidate) & isfinite(obj_candidate);
            if ~isempty(fidelity_file)
                save(fidelity_file, 'fidelity', '-v7.3');
            end
            Out_temp_dec = [Out_dec; solution_candidate_filtered]; 
            Out_temp_obj = [Out_obj; obj_candidate];
            
%             mask = any(Out_temp_dec ~= 0, 2);  % 标记非全0的个体
%             Out_temp_dec = Out_temp_dec(mask, :);
%             Out_temp_obj = Out_temp_obj(mask, :);
            
            TEMP = [Out_temp_dec, Out_temp_obj];
            FEsUsed = FEsUsed + size(solution_candidate_filtered,1);
        else
            Out_temp_dec = Out_dec;
            Out_temp_obj = Out_obj;
            
%             mask = any(Out_temp_dec ~= 0, 2);  % 标记非全0的个体
%             Out_temp_dec = Out_temp_dec(mask, :);
%             Out_temp_obj = Out_temp_obj(mask, :);
            
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
%        %% bio求解之前需要删掉全为0的个体否则影响个体求解
%         nonzero_mask = any(archive_best_dec ~= 0, 2);  % 标记非全0的个体
%         archive_best_dec = archive_best_dec(nonzero_mask, :);
%         archive_best_obj = archive_best_obj(nonzero_mask, :);
        %%  FE加速比求解     
        if isempty(reachedFE)
            if isKey(threshold_map, name)
                target_threshold = threshold_map(name);
            else
                error('未定义 %s 的阈值', name);
            end
            if any(all(archive_best_obj <= target_threshold, 2))
                reachedFE = FEsUsed;  % 记录首次满足条件的 FE
            end
        end

%        %% ===========测试函数计算指标============
%         IGD = calIGD(archive_best_obj, optimum);
%         if IGD < best_IGD
%             best_IGD = IGD;
%         end
%         % 记录IGD历史
%         igd_history = [igd_history, IGD];
%         % 计算滑动窗口内IGD变化均值
%         if length(igd_history) >= window_size + 1
%             recent_changes = abs(diff(igd_history(end-window_size:end)));
%             avg_change = mean(recent_changes);
%         else
%             avg_change = inf; % 不足窗口时不触发切换
%         end
%         switch phase_mode
%             case 'self'
%                 if avg_change < improvement_threshold && (FEsUsed - last_switch_FE) >= cooldown_FE
%                     fprintf('>>> IGD stable, switch to transfer phase.\n');
%                     phase_mode = 'transfer'; % transfer
%                     transfer_attempts = 0;
%                     last_transfer_IGD = best_IGD;
%                     last_switch_FE = FEsUsed;  % 记录这次切换时的FEs
%                 end
%               case 'transfer'
%                  if IGD < last_transfer_IGD - 1e-6  % 有改进
%                     fprintf('>>> Transfer improving, stay in transfer phase.\n');
% %                     phase_mode = 'self';
%                     transfer_attempts = 0;
%                     last_transfer_IGD = IGD;
%                 else
%                     transfer_attempts = transfer_attempts + 1;
%                     fprintf('>>> Transfer no improvement (%d/%d).\n', transfer_attempts, max_transfer_attempts);
%                     if transfer_attempts >= max_transfer_attempts
%                         fprintf('>>> Transfer abandoned. Switch to self.\n');
%                         phase_mode = 'self';
%                         last_switch_FE = FEsUsed;  % 记录切回self的时间
%                     end
%                  end
%         end
%         IGD_record = [IGD_record; best_IGD];
%         FE_record = [FE_record; FEsUsed];
% %         if paras.ada_flag == 0
%         fprintf([name_saea,'-CKT-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,best_IGD);
        
        %% ===========bio应用计算指标============
        % 初始化 min_vals 和 max_vals（可在外层初始化）
        if ~norm_flag
            min_vals = min(archive_best_obj, [], 1);
            max_vals = max(archive_best_obj, [], 1);
            norm_flag = true;  % 标记已经初始化过
        end
        % 动态更新归一化上下界
        min_vals = min(min_vals, min(archive_best_obj, [], 1));
        max_vals = max(max_vals, max(archive_best_obj, [], 1));
        range_vals = max(max_vals - min_vals, 1e-8);  % 避免除以 0
        % 对整个 archive_best_obj 中的所有解做归一化
        norm_objs = (archive_best_obj - min_vals) ./ range_vals;  % size: NxM
        % 计算每个解的平均归一化目标值
        avg_per_solution = mean(norm_objs, 2);  % size: Nx1
        % 取最小的平均值作为当前代表指标（你也可以改为 mean(...)）
        avg_norm_obj = mean(avg_per_solution);  %
        all_avg_norm = [all_avg_norm, avg_norm_obj];
        % 记录指标变化历史
        norm_obj_history = [norm_obj_history, avg_norm_obj];
        % 滑动窗口判断是否稳定
        if length(norm_obj_history) >= window_size + 1
            recent_changes = abs(diff(norm_obj_history(end - window_size : end)));
            avg_change = mean(recent_changes);
        else
            avg_change = inf;  % 初期不触发
        end
        % ===========阶段切换逻辑============
        switch phase_mode
            case 'self' %  && (FEsUsed - last_switch_FE) >= cooldown_FE
                if avg_change < improvement_threshold
                    fprintf('>>> Normalized Obj stable, switch to transfer phase.\n');
                    phase_mode = 'transfer';
                    transfer_attempts = 0;
                    last_transfer_obj = avg_norm_obj;
                    last_switch_FE = FEsUsed;
                end
            case 'transfer'
                if avg_norm_obj < last_transfer_obj - 1e-4
                    fprintf('>>> Transfer improving, stay in transfer phase.\n');
                    transfer_attempts = 0;
                    last_transfer_obj = avg_norm_obj;
                else
                    transfer_attempts = transfer_attempts + 1;
                    fprintf('>>> Transfer no improvement (%d/%d).\n', transfer_attempts, max_transfer_attempts);
                    if transfer_attempts >= max_transfer_attempts
                        fprintf('>>> Transfer abandoned. Switch to self.\n');
                        phase_mode = 'self';
                        last_switch_FE = FEsUsed;
                    end
                end
        end
        % 保存记录
        FE_record = [FE_record; FEsUsed];
        norm_obj_record = [norm_obj_record; avg_norm_obj];
        % 打印当前状态
        fprintf([name_saea, '-CKT-FEsUsed: %d, norm_avg_obj: %.4f\n'], FEsUsed, avg_norm_obj);
%         else
%             fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,IGD);
%         end
    end
    runtime = toc;
    
%     %% 输出保存的模型结果 === RMSE: 模型 A vs B 对比 (针对目标 m) ===
%     save('model_compare_results.mat', ...
%     'rmseA_record', 'rmseB_record', ...
%     'maeA_record', 'maeB_record', ...
%     'R2A_record', 'R2B_record', ...
%     'FE_Record');
%     M = size(rmseA_record, 2);
%     for m = 1:M
%         figure;
%         plot(FE_record, rmseA_record(:, m), 'b*--', 'LineWidth', 1.8); hold on;
%         plot(FE_record, rmseB_record(:, m), 'r*-', 'LineWidth', 1.8);
%         xlabel('Function Evaluations', 'FontSize', 12, 'FontName', 'Times New Roman');
%         ylabel(sprintf('RMSE (Objective %d)', m), 'FontSize', 12, 'FontName', 'Times New Roman');
%         title(sprintf('RMSE Comparison on Objective %d', m), 'FontSize', 13, 'FontName', 'Times New Roman');
%         legend('Original Model', 'Warm-start Model', 'FontSize', 11, 'FontName', 'Times New Roman', 'Location', 'best');
%         xlim([0, 1000]);
%         grid on;
%         set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
%         set(gcf, 'Units', 'inches', 'Position', [1, 1, 5.2, 4]);
%         % 保存 RMSE 图
%         print(gcf, sprintf('results/RMSE_Obj%d_WFG1.png', m), '-dpng', '-r600');
%         % === R2: 模型 A vs B 对比 (针对目标 m) ===
%         figure;
%         plot(FE_record, R2A_record(:, m), 'b*--', 'LineWidth', 1.8); hold on;
%         plot(FE_record, R2B_record(:, m), 'r*-', 'LineWidth', 1.8);
%         xlabel('Function Evaluations', 'FontSize', 12, 'FontName', 'Times New Roman');
%         ylabel(sprintf('R^{2}(Objective %d)', m), 'FontSize', 12, 'FontName', 'Times New Roman');
%         title(sprintf('R^{2} Comparison on Objective %d', m), 'FontSize', 13, 'FontName', 'Times New Roman');
%         legend('Original Model', 'Warm-start Model', 'FontSize', 11, 'FontName', 'Times New Roman', 'Location', 'best');
%         xlim([0, 1000]);
%         grid on;
%         set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
%         set(gcf, 'Units', 'inches', 'Position', [1, 1, 5.2, 4]);
%         % 保存 MAE 图
%         print(gcf, sprintf('results/R2_Obj%d_WFG1.png', m), '-dpng', '-r600');
%     end

    
    %% =====输出archive=======
    Out_dec = archive_best_dec;
    Out_obj = archive_best_obj;
    % ===== 拥挤度截断控制 archive大小 =====
%     archive_limit = 500;
    N = size(Out_obj, 1);
    if N > archive_limit
        dis = CrowdingDistance(Out_obj);  % 计算拥挤距离
        [~, idx] = sort(dis, 'descend');  % 根据拥挤距离从大到小排序
        idx = idx(1:archive_limit);       % 保留前 archive_limit 个解
        Out_dec = Out_dec(idx, :);
        Out_obj = Out_obj(idx, :);
    end
    fidelity.optimizer = name_saea;
    fidelity.FE_final = FEsUsed;
    fidelity.iterations = iteration;
    fidelity.runtime = runtime;
    fidelity.FE_record = FE_record;
    fidelity.norm_obj_record = norm_obj_record;
    fidelity.completed = true;
    if ~isempty(fidelity_file)
        save(fidelity_file, 'fidelity', '-v7.3');
    end
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

function dom_flag = is_dominated(a, B)
% 判断解 a 是否被 B 中任一解支配
    dom_flag = any(all(B <= a, 2) & any(B < a, 2));
end

function idx_keep = filter_nondominated(A, B)
% 保留 A 中未被 B 支配的个体索引
    idx_keep = [];
    for i = 1:size(A,1)
        if ~is_dominated(A(i,:), B)
            idx_keep = [idx_keep; i];
        end
    end
end

% function mask = epsilon_dominance_relation(popObj, refObj, epsilon)
% % EPSILON_DOMINANCE_RELATION 判定 popObj 中哪些解 ε-支配 refObj 中的任意解
% %
% % 输入:
% %   popObj   - N x M 个体目标函数值矩阵
% %   refObj   - K x M 参考解目标函数值矩阵（如 archive_best_obj）
% %   epsilon  - ε 支配的松弛因子（标量）
% %
% % 输出:
% %   mask     - 逻辑向量，表示 popObj 中哪些个体 ε-支配 refObj 中的任一解
% 
%     N = size(popObj, 1);
%     K = size(refObj, 1);
%     mask = false(N, 1);  % 初始化掩码
% 
%     for i = 1:N
%         % 该解是否支配 refObj 中的任意一个解（存在）
%         for j = 1:K
%             % 判断是否满足 ε-支配条件
%             if all(popObj(i,:) <= refObj(j,:) + epsilon) && ...
%                any(popObj(i,:) < refObj(j,:) - epsilon)
%                 mask(i) = true;
%                 break;  % 当前个体已支配任一解，跳出内层循环
%             end
%         end
%     end
% end




%% =================原版===============
% function [Out_dec,Out_obj, metric, all_avg_norm, reachedFE,runtime]= sas_ckt(fun,lb,ub,FEsMax,name_saea,knowledge_base,paras,M,optimum,name)
% 
%     % =================Phase 1: Initialization===================
%     if isempty(name_saea)
%         name_saea = 'NSGAII';
%     end
% %     num_sources = length(knowledge_base); 
% %     interval_transfer = paras.interval;
%     num_initial_solutions = size(paras.pop,1); % 初始化种群个数
%     solutions_initial = paras.pop;   % 初始种群
%     objs_initial = zeros(num_initial_solutions,M);
%     %% 如果后面变量需要迁移的话，这里就不需要计算目标函数值
%     for i = 1:num_initial_solutions
%         objs_initial(i,:) = fun(solutions_initial(i,:));  %fun是传递的目标任务
%     end
%     FEsUsed = num_initial_solutions;
%     database = [solutions_initial, objs_initial]; % 有目标值的种群
% %     dim = length(lb);
%     transfer_states = [];
%     if strcmp(name_saea,'AutoSAEA')
%         load('optimizers\auto-saea\autosaea_paras.mat');
%         options_saea.func_exp = fun;
%     else
%         options_saea = [];
%     end
% 
%     Out_dec = database(:,1:end-M); 
%     Out_obj = database(:,end-M+1:end);
% 
%     %% 设置一个存档保存非支配解作为最后的输出
%     [~,FrontNo] = NSGASort(Out_obj);
%     index = FrontNo==1;
%     archive_Dec = Out_dec(index,:);
%     archive_obj = Out_obj(index,:);
% 
%     IGD = calIGD(archive_obj, optimum);
%     fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,IGD);
% 
%     %% 初始化存档参数设置
%     best_IGD = inf;
%     phase_mode = 'self'; 
%     igd_last = inf;  % 初始化前一代 IGD
%     isFirstWarmStart = true;  % 标志是否第一次进入
% 
%     norm_flag = false;
%     norm_obj_history = [];
%     norm_obj_record = [];
%     all_avg_norm = [];
% 
%     %% 迁移转换初始化    
%     igd_history = [];
%     window_size = 10;
%     improvement_threshold = 1e-4;
%     max_transfer_attempts = 5;
%     transfer_attempts = 0;
%     cooldown_FE = 50;        % 冷却阶段，单位是函数评估次数（可调）
%     last_switch_FE = 0;       % 上次阶段切换的 FEs
% 
%     IGD_record = [];
%     FE_record = []; 
%     archive_best_obj = archive_obj;
%     archive_best_dec = archive_Dec;
%     archive_limit = 500; 
% 
%     % 加速比求解
%     reachedFE = [];
%     threshold_map = containers.Map( ...
%     {'bio1', 'bio2', 'bio3', 'bio4','bio5','bio6','bio7','bio8','bio9','bio10'}, ...
%     {[0.3, 3], [0.3, 4], [0.32, 4], [0.29, 10], [0.36, 7], [0.34, 5], [0.35, 2], [0.3, 25], [0.3, 10], [0.3, 11]} ...
% );
% 
% %     %% 记录模型的性能提升初始化集合
% %     rmseA_record = [];
% %     rmseB_record = [];
% %     maeA_record = [];
% %     maeB_record = [];
% %     FE_Record = [];
% %     R2A_record = []; 
% %     R2B_record = [];
% 
%     % =========================Phase 2: Optimization=====================
% %     %% 测试集生成验证模型性能
% %     Ytest = zeros(num_initial_solutions,M);
% %     num_samples = 100;      % 测试集样本数量
% %     num_vars = 20;          % 决策变量维度
% %     lower_bounds = zeros(1, num_vars); 
% %     upper_bounds = ones(1, num_vars);   
% %     % Latin Hypercube Sampling
% %     lhs_unit = lhsdesign(num_samples, num_vars);  % 在 [0,1]^n 中采样
% %     X_test = repmat(lower_bounds, num_samples, 1) + ...
% %      lhs_unit .* repmat((upper_bounds - lower_bounds), num_samples, 1);
% %      for i = 1:num_samples
% %         Ytest(i,:) = fun(X_test(i,:)); 
% %      end
%     tic;
%     while FEsUsed < FEsMax
% %          % ------未迁移之前的数据集------
% %         target_pop = database(:,1:end-M);
% %         objs = database(:,end-M+1:end);
% %         func_surrogateA.gpr = cell(1, M);
% %         for i = 1:M
% %             func_surrogateA.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
% %         end
% %         for i = 1:M
% %             [objs_pre(:, i),objs_std(:, i)] = func_surrogateA.gpr{i}.func(X_test);
% %         end
% %         Y_predA = objs_pre;
%        % Step 1: 训练/更新代理模型 
%          if isFirstWarmStart
%             [newdatabase, FEs_inc, ~, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
%             % [newdatabase, FEs_inc, ~, transfer_source] = warm_start_variable_transfer_LAE(database, knowledge_base,fun, M);
% 
% 
%             database = newdatabase;
%             FEsUsed = FEsUsed + FEs_inc;
%             isFirstWarmStart = false;
%          end    
%         target_pop = database(:,1:end-M);
%         objs = database(:,end-M+1:end);
%         func_surrogate.gpr = cell(1, M);
%         for i = 1:M
%             func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
%         end
%         if ~isFirstWarmStart
%             % Step 2: 评估代理模型误差
%             current_rmse = evaluate_model_error(func_surrogate, database,M);
%             disp(current_rmse);
%             %% Step 3: 判断是否允许迁移 代理模型的精度这里需要消融
% %             threshold_rmse = 0.05;  
%             threshold_rmse = 1e-4*5; 
%             if current_rmse > threshold_rmse
%                 % 代理模型还不够准，继续warm start
%                 [newdatabase, FEs_inc, success, transfer_source] = warm_start_variable_transfer(database, knowledge_base,fun, M);
%                 database = newdatabase;
%                 FEsUsed = FEsUsed + FEs_inc;
%                 if success
%                     update_model_flag = true;
%                     % warm start 后更新 surrogate 模型
%                     target_pop = database(:,1:end-M);
%                     objs = database(:,end-M+1:end);
%                     for i = 1:M
%                         func_surrogate.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
%                     end
%                 end
%             end
%         end
% %         %% =====计算 RMSE R2 MAE=====
% %         target_pop = database(:,1:end-M);
% %         objs = database(:,end-M+1:end);
% %         func_surrogateB.gpr = cell(1, M);
% %         for i = 1:M
% %             func_surrogateB.gpr{i} = surrogate_model(target_pop, objs(:, i), 'gpr', M);
% %         end
% %         for i = 1:M
% %             [objs_pre(:, i),objs_std(:, i)] = func_surrogateB.gpr{i}.func(X_test);
% %         end
% %         Y_predB = objs_pre;
% %         rmseA = sqrt(mean((Y_predA - Ytest).^2));
% %         rmseB = sqrt(mean((Y_predB - Ytest).^2));
% %         maeA = mean(abs(Y_predA - Ytest));
% %         maeB = mean(abs(Y_predB - Ytest));
% %         % R2 for Model A
% %         SS_res_A = sum((Ytest - Y_predA).^2);                  % Residual sum of squares
% %         SS_tot_A = sum((Ytest - mean(Ytest)).^2);              % Total sum of squares
% %         R2_A = 1 - SS_res_A ./ SS_tot_A;                       % R^2 per objective
% %         % R2 for Model B
% %         SS_res_B = sum((Ytest - Y_predB).^2);
% %         SS_tot_B = sum((Ytest - mean(Ytest)).^2);
% %         R2_B = 1 - SS_res_B ./ SS_tot_B;
% %         % 记录
% %         rmseA_record(end+1,:) = rmseA;
% %         rmseB_record(end+1,:) = rmseB;
% %         maeA_record(end+1,:) = maeA;
% %         maeB_record(end+1,:) = maeB;
% %         R2A_record(end+1, :) = R2_A;
% %         R2B_record(end+1, :) = R2_B;
% %         FE_Record(end+1) = FEsUsed;
% %         % 打印
% %         fprintf('FEs = %d | RMSE A: %.4f, B: %.4f | MAE A: %.4f, B: %.4f\n', ...
% %             FEsUsed, mean(rmseA), mean(rmseB), mean(maeA), mean(maeB));
% 
%         %% 自迁移
%         [solution_in,improvement_in,options_saea] = self_opt(database,name_saea,options_saea,func_surrogate,M,archive_best_obj);
% %         solution_candidate = solution_in;
% %         solution_candidate_obj = improvement_in;
%         %% 动态迁移强度设置
%         if strcmp(phase_mode, 'transfer')
%             [solution_ex,improvement_ex] = transfer_opt(database,name_saea,options_saea,func_surrogate,M,transfer_source,knowledge_base);
%             idx_B = filter_nondominated(improvement_ex, improvement_in);
%             solution_candidate = solution_ex(idx_B, :);
%             solution_candidate_obj = improvement_ex(idx_B, :);
%         else
%             solution_candidate = solution_in;
%             solution_candidate_obj = improvement_in;
%         end
%         solution_candidate_filtered = solution_candidate;
% 
%        %% 填充个体数
%         K = 11;
%         num = min(K, size(solution_candidate_filtered,1));
%         selected_idx = randperm(size(solution_candidate_filtered,1), num);
%         solution_candidate_filtered = solution_candidate_filtered(selected_idx, :);
%         if ~isempty(solution_candidate_filtered)
%             obj_candidate = fun(solution_candidate_filtered);
%             Out_temp_dec = [Out_dec; solution_candidate_filtered]; 
%             Out_temp_obj = [Out_obj; obj_candidate];
% 
% %             mask = any(Out_temp_dec ~= 0, 2);  % 标记非全0的个体
% %             Out_temp_dec = Out_temp_dec(mask, :);
% %             Out_temp_obj = Out_temp_obj(mask, :);
% 
%             TEMP = [Out_temp_dec, Out_temp_obj];
%             FEsUsed = FEsUsed + size(solution_candidate_filtered,1);
%         else
%             Out_temp_dec = Out_dec;
%             Out_temp_obj = Out_obj;
% 
% %             mask = any(Out_temp_dec ~= 0, 2);  % 标记非全0的个体
% %             Out_temp_dec = Out_temp_dec(mask, :);
% %             Out_temp_obj = Out_temp_obj(mask, :);
% 
%             TEMP = [Out_temp_dec, Out_temp_obj];
%         end
% 
%         % =======种群更新======
%         [rank,~,~] = NSGASort(Out_temp_obj);
%         database = TEMP(rank(1:num_initial_solutions),:);
%         % ======更新archive=archive+population======
%         combined_dec = [archive_best_dec; Out_temp_dec];
%         combined_obj = [archive_best_obj; Out_temp_obj];
%         % 去除支配解，保留真正的非支配个体
%         [~, FrontNo_combined] = NSGASort(combined_obj);
%         idx_best = FrontNo_combined == 1;
%         % 提取合并后的非支配解
%         archive_best_dec = combined_dec(idx_best, :);
%         archive_best_obj = combined_obj(idx_best, :);
% %        %% bio求解之前需要删掉全为0的个体否则影响个体求解
% %         nonzero_mask = any(archive_best_dec ~= 0, 2);  % 标记非全0的个体
% %         archive_best_dec = archive_best_dec(nonzero_mask, :);
% %         archive_best_obj = archive_best_obj(nonzero_mask, :);
%         %%  FE加速比求解     
%         if isempty(reachedFE)
%             if isKey(threshold_map, name)
%                 target_threshold = threshold_map(name);
%             else
%                 error('未定义 %s 的阈值', name);
%             end
%             if any(all(archive_best_obj <= target_threshold, 2))
%                 reachedFE = FEsUsed;  % 记录首次满足条件的 FE
%             end
%         end
% 
% %        %% ===========测试函数计算指标============
% %         IGD = calIGD(archive_best_obj, optimum);
% %         if IGD < best_IGD
% %             best_IGD = IGD;
% %         end
% %         % 记录IGD历史
% %         igd_history = [igd_history, IGD];
% %         % 计算滑动窗口内IGD变化均值
% %         if length(igd_history) >= window_size + 1
% %             recent_changes = abs(diff(igd_history(end-window_size:end)));
% %             avg_change = mean(recent_changes);
% %         else
% %             avg_change = inf; % 不足窗口时不触发切换
% %         end
% %         switch phase_mode
% %             case 'self'
% %                 if avg_change < improvement_threshold && (FEsUsed - last_switch_FE) >= cooldown_FE
% %                     fprintf('>>> IGD stable, switch to transfer phase.\n');
% %                     phase_mode = 'transfer'; % transfer
% %                     transfer_attempts = 0;
% %                     last_transfer_IGD = best_IGD;
% %                     last_switch_FE = FEsUsed;  % 记录这次切换时的FEs
% %                 end
% %               case 'transfer'
% %                  if IGD < last_transfer_IGD - 1e-6  % 有改进
% %                     fprintf('>>> Transfer improving, stay in transfer phase.\n');
% % %                     phase_mode = 'self';
% %                     transfer_attempts = 0;
% %                     last_transfer_IGD = IGD;
% %                 else
% %                     transfer_attempts = transfer_attempts + 1;
% %                     fprintf('>>> Transfer no improvement (%d/%d).\n', transfer_attempts, max_transfer_attempts);
% %                     if transfer_attempts >= max_transfer_attempts
% %                         fprintf('>>> Transfer abandoned. Switch to self.\n');
% %                         phase_mode = 'self';
% %                         last_switch_FE = FEsUsed;  % 记录切回self的时间
% %                     end
% %                  end
% %         end
% %         IGD_record = [IGD_record; best_IGD];
% %         FE_record = [FE_record; FEsUsed];
% % %         if paras.ada_flag == 0
% %         fprintf([name_saea,'-CKT-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,best_IGD);
% 
%         %% ===========bio应用计算指标============
%         % 初始化 min_vals 和 max_vals（可在外层初始化）
%         if ~norm_flag
%             min_vals = min(archive_best_obj, [], 1);
%             max_vals = max(archive_best_obj, [], 1);
%             norm_flag = true;  % 标记已经初始化过
%         end
%         % 动态更新归一化上下界
%         min_vals = min(min_vals, min(archive_best_obj, [], 1));
%         max_vals = max(max_vals, max(archive_best_obj, [], 1));
%         range_vals = max(max_vals - min_vals, 1e-8);  % 避免除以 0
%         % 对整个 archive_best_obj 中的所有解做归一化
%         norm_objs = (archive_best_obj - min_vals) ./ range_vals;  % size: NxM
%         % 计算每个解的平均归一化目标值
%         avg_per_solution = mean(norm_objs, 2);  % size: Nx1
%         % 取最小的平均值作为当前代表指标（你也可以改为 mean(...)）
%         avg_norm_obj = mean(avg_per_solution);  %
%         all_avg_norm = [all_avg_norm, avg_norm_obj];
%         % 记录指标变化历史
%         norm_obj_history = [norm_obj_history, avg_norm_obj];
%         % 滑动窗口判断是否稳定
%         if length(norm_obj_history) >= window_size + 1
%             recent_changes = abs(diff(norm_obj_history(end - window_size : end)));
%             avg_change = mean(recent_changes);
%         else
%             avg_change = inf;  % 初期不触发
%         end
%         % ===========阶段切换逻辑============
%         switch phase_mode
%             case 'self' %  && (FEsUsed - last_switch_FE) >= cooldown_FE
%                 if avg_change < improvement_threshold
%                     fprintf('>>> Normalized Obj stable, switch to transfer phase.\n');
%                     phase_mode = 'transfer';
%                     transfer_attempts = 0;
%                     last_transfer_obj = avg_norm_obj;
%                     last_switch_FE = FEsUsed;
%                 end
%             case 'transfer'
%                 if avg_norm_obj < last_transfer_obj - 1e-4
%                     fprintf('>>> Transfer improving, stay in transfer phase.\n');
%                     transfer_attempts = 0;
%                     last_transfer_obj = avg_norm_obj;
%                 else
%                     transfer_attempts = transfer_attempts + 1;
%                     fprintf('>>> Transfer no improvement (%d/%d).\n', transfer_attempts, max_transfer_attempts);
%                     if transfer_attempts >= max_transfer_attempts
%                         fprintf('>>> Transfer abandoned. Switch to self.\n');
%                         phase_mode = 'self';
%                         last_switch_FE = FEsUsed;
%                     end
%                 end
%         end
%         % 保存记录
%         FE_record = [FE_record; FEsUsed];
%         norm_obj_record = [norm_obj_record; avg_norm_obj];
%         % 打印当前状态
%         fprintf([name_saea, '-CKT-FEsUsed: %d, norm_avg_obj: %.4f\n'], FEsUsed, avg_norm_obj);
% %         else
% %             fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, best_igd: %.2f\n'],FEsUsed,IGD);
% %         end
%     end
%     runtime = toc;
% 
% %     %% 输出保存的模型结果 === RMSE: 模型 A vs B 对比 (针对目标 m) ===
% %     save('model_compare_results.mat', ...
% %     'rmseA_record', 'rmseB_record', ...
% %     'maeA_record', 'maeB_record', ...
% %     'R2A_record', 'R2B_record', ...
% %     'FE_Record');
% %     M = size(rmseA_record, 2);
% %     for m = 1:M
% %         figure;
% %         plot(FE_record, rmseA_record(:, m), 'b*--', 'LineWidth', 1.8); hold on;
% %         plot(FE_record, rmseB_record(:, m), 'r*-', 'LineWidth', 1.8);
% %         xlabel('Function Evaluations', 'FontSize', 12, 'FontName', 'Times New Roman');
% %         ylabel(sprintf('RMSE (Objective %d)', m), 'FontSize', 12, 'FontName', 'Times New Roman');
% %         title(sprintf('RMSE Comparison on Objective %d', m), 'FontSize', 13, 'FontName', 'Times New Roman');
% %         legend('Original Model', 'Warm-start Model', 'FontSize', 11, 'FontName', 'Times New Roman', 'Location', 'best');
% %         xlim([0, 1000]);
% %         grid on;
% %         set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
% %         set(gcf, 'Units', 'inches', 'Position', [1, 1, 5.2, 4]);
% %         % 保存 RMSE 图
% %         print(gcf, sprintf('results/RMSE_Obj%d_WFG1.png', m), '-dpng', '-r600');
% %         % === R2: 模型 A vs B 对比 (针对目标 m) ===
% %         figure;
% %         plot(FE_record, R2A_record(:, m), 'b*--', 'LineWidth', 1.8); hold on;
% %         plot(FE_record, R2B_record(:, m), 'r*-', 'LineWidth', 1.8);
% %         xlabel('Function Evaluations', 'FontSize', 12, 'FontName', 'Times New Roman');
% %         ylabel(sprintf('R^{2}(Objective %d)', m), 'FontSize', 12, 'FontName', 'Times New Roman');
% %         title(sprintf('R^{2} Comparison on Objective %d', m), 'FontSize', 13, 'FontName', 'Times New Roman');
% %         legend('Original Model', 'Warm-start Model', 'FontSize', 11, 'FontName', 'Times New Roman', 'Location', 'best');
% %         xlim([0, 1000]);
% %         grid on;
% %         set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
% %         set(gcf, 'Units', 'inches', 'Position', [1, 1, 5.2, 4]);
% %         % 保存 MAE 图
% %         print(gcf, sprintf('results/R2_Obj%d_WFG1.png', m), '-dpng', '-r600');
% %     end
% 
% 
%     %% =====输出archive=======
%     Out_dec = archive_best_dec;
%     Out_obj = archive_best_obj;
%     % ===== 拥挤度截断控制 archive大小 =====
% %     archive_limit = 500;
%     N = size(Out_obj, 1);
%     if N > archive_limit
%         dis = CrowdingDistance(Out_obj);  % 计算拥挤距离
%         [~, idx] = sort(dis, 'descend');  % 根据拥挤距离从大到小排序
%         idx = idx(1:archive_limit);       % 保留前 archive_limit 个解
%         Out_dec = Out_dec(idx, :);
%         Out_obj = Out_obj(idx, :);
%     end
%     metric = [FE_record,IGD_record];
% end
% 
% function dis = CrowdingDistance(obj)
%     [N, M] = size(obj);
%     dis = zeros(N,1);
%     Fmax = max(obj,[],1);
%     Fmin = min(obj,[],1);
%     for i = 1:M
%         [~, rank] = sort(obj(:,i));
%         dis(rank(1)) = inf;
%         dis(rank(end)) = inf;
%         for j = 2:N-1
%             if Fmax(i) - Fmin(i) == 0
%                 d = 0;
%             else
%                 d = (obj(rank(j+1),i) - obj(rank(j-1),i)) / (Fmax(i) - Fmin(i));
%             end
%             dis(rank(j)) = dis(rank(j)) + d;
%         end
%     end
% end
% 
% function dom_flag = is_dominated(a, B)
% % 判断解 a 是否被 B 中任一解支配
%     dom_flag = any(all(B <= a, 2) & any(B < a, 2));
% end
% 
% function idx_keep = filter_nondominated(A, B)
% % 保留 A 中未被 B 支配的个体索引
%     idx_keep = [];
%     for i = 1:size(A,1)
%         if ~is_dominated(A(i,:), B)
%             idx_keep = [idx_keep; i];
%         end
%     end
% end
% 
% % function mask = epsilon_dominance_relation(popObj, refObj, epsilon)
% % % EPSILON_DOMINANCE_RELATION 判定 popObj 中哪些解 ε-支配 refObj 中的任意解
% % %
% % % 输入:
% % %   popObj   - N x M 个体目标函数值矩阵
% % %   refObj   - K x M 参考解目标函数值矩阵（如 archive_best_obj）
% % %   epsilon  - ε 支配的松弛因子（标量）
% % %
% % % 输出:
% % %   mask     - 逻辑向量，表示 popObj 中哪些个体 ε-支配 refObj 中的任一解
% % 
% %     N = size(popObj, 1);
% %     K = size(refObj, 1);
% %     mask = false(N, 1);  % 初始化掩码
% % 
% %     for i = 1:N
% %         % 该解是否支配 refObj 中的任意一个解（存在）
% %         for j = 1:K
% %             % 判断是否满足 ε-支配条件
% %             if all(popObj(i,:) <= refObj(j,:) + epsilon) && ...
% %                any(popObj(i,:) < refObj(j,:) - epsilon)
% %                 mask(i) = true;
% %                 break;  % 当前个体已支配任一解，跳出内层循环
% %             end
% %         end
% %     end
% % end
% 
