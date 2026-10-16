% function [database, FEsUsed_add, success, source] = warm_start_variable_transfer(database, knowledge_base, fun, M)
%     success = false;
%     FEsUsed_add = 0;
%     try
%         target_pop = database(:, 1:end-M);
%         target_obj = database(:, end-M+1:end);
%         dimt = size(target_pop, 2);
%         source_pops = {knowledge_base.solutions};
%         mind = cellfun(@(x) size(x, 2), source_pops);
%         dims = min(mind);
%         proj_dim = min(dimt, dims);
% 
%         % Step 1: Find most similar source task 这里是随机投影+MMD
%         most_similar_idx_history = [];
%         [most_similar_idx, ~] = FindMostSimilarTask(target_pop, source_pops, proj_dim);
%         most_similar_idx_history(end+1) = most_similar_idx;   % 记录这次的索引
%         disp(most_similar_idx_history);
%         transfer_source = knowledge_base(most_similar_idx);
%         idx = transfer_source.key_var;  %var_d
%         transfer_var = mean(transfer_source.solutions(:, idx)); % 这里是取了均值的
% 
%         % Step 2: 从最差的帕累托前沿开始，逐层选择个体，直到达到替换数量
%         pop_size = size(target_pop, 1);
%         [~, FrontNo] = NSGASort(target_obj);
%         max_front = max(FrontNo);
%         %% =======这里有一个消融变量，即选择变量替换的比列======
%         num_replace = ceil(0.2 * pop_size);
%         selected_ind = [];  % 初始化选择个体下标
%         for front = max_front:-1:1
%             inds = find(FrontNo == front);
%             selected_ind = [selected_ind; inds(:)];
%             if length(selected_ind) >= num_replace
%                 selected_ind = selected_ind(1:num_replace);
%                 break;
%             end
%         end
%         if length(selected_ind) < num_replace
%             remain = setdiff(1:pop_size, selected_ind);
%             extra = remain(randperm(length(remain), num_replace - length(selected_ind)));
%             selected_ind = [selected_ind; extra(:)];
%         end
% 
%         var = transfer_var;
%         selected_idx = idx;
%         % Step 2.5: 去除超过目标维度的变量
%         valid_mask = selected_idx <= dimt;
%         selected_idx = selected_idx(valid_mask);
%         var = var(valid_mask);
%         % 如果源任务的变量数小于目标任务的维度
%         if length(var) < dimt
%             % 不需要做任何改变，直接用源任务的变量
%             selected_idx = idx(1:length(var));  % 根据源任务变量的数量限制索引
%         end
%         if length(var) > dimt
%             % 如果源任务变量数大于目标维度，随机选择 dimt 个索引
%             selected_idx = idx(randperm(length(idx), dimt));
%             var = transfer_var(selected_idx);  % 确保替换的变量维度一致
%         end
% 
%         % Step 5: 替换目标种群中的个体
%         for i = 1:length(selected_ind)
%             ind = selected_ind(i);
%             % 检查 selected_idx 维度和 var 维度是否匹配
%             if length(selected_idx) ~= length(var)
%                 error('Selected index dimensions do not match transfer variable length.');
%             end
%             % 进行替换，只替换目标种群中指定维度的变量
%             target_pop(ind, selected_idx) = var;
%             target_obj(ind, :) = fun(target_pop(ind, :));  % 重新计算目标值
%         end
% 
%         % Step 4: Update database with modified individuals
%         database = [target_pop, target_obj];
%         FEsUsed_add = length(selected_ind); 
%         success = true;
%         source = transfer_source;
%     catch
%         % 异常时输出参数设为默认值
%         database = [target_pop, target_obj];
%         FEsUsed_add = 0;
%         success = false;
%         source = [];  
%         return;
%     end
% end

function [database, FEsUsed_add, success, source] = ...
    warm_start_variable_transfer(database, knowledge_base, fun, M)

    success = false;
    FEsUsed_add = 0;
    source = [];

    if isempty(knowledge_base)
        warning('STOK:EmptyKB', '候选源任务知识库为空，跳过迁移。');
        return;
    end

    target_pop = database(:, 1:end-M);
    target_obj = database(:, end-M+1:end);
    dimt = size(target_pop, 2);
    pop_size = size(target_pop, 1);

    if pop_size == 0
        return;
    end

    % Step 1: 使用原来的随机投影 + MMD 选择最相似源任务
    % knowledge_base 由 main 传入，已排除当前目标任务。
    source_pops = {knowledge_base.solutions};
    source_dims = cellfun(@(x) size(x, 2), source_pops);
    proj_dim = min(dimt, min(source_dims));

    [most_similar_idx, ~] = ...
        FindMostSimilarTask(target_pop, source_pops, proj_dim);

    assert(isnumeric(most_similar_idx) && ...
           isscalar(most_similar_idx) && ...
           isfinite(most_similar_idx) && ...
           most_similar_idx == fix(most_similar_idx) && ...
           most_similar_idx >= 1 && ...
           most_similar_idx <= numel(knowledge_base), ...
           'STOK:SourceIndex', 'MMD 返回的源任务下标无效。');

    transfer_source = knowledge_base(most_similar_idx);

    % 此处为候选库局部下标，不是真实 Bio 编号
    fprintf('Full STOK: selected source local index = %d\n', ...
        most_similar_idx);

    % 保留选中源任务的关键变量选择机制
    idx = transfer_source.key_var;
    idx = idx(:)';

    % 同步过滤索引，避免越界
    source_dim = size(transfer_source.solutions, 2);
    valid_mask = isfinite(idx) & idx >= 1 & ...
                 idx == fix(idx) & idx <= dimt & idx <= source_dim;
    selected_idx = idx(valid_mask);

    if isempty(selected_idx) || isempty(transfer_source.solutions)
        warning('STOK:NoTransferVariables', ...
            '选中源任务没有可迁移的有效变量，跳过本次迁移。');
        return;
    end

    % 保留原来的均值迁移方式
    transfer_var = mean(transfer_source.solutions(:, selected_idx), 1);

    % Step 2: 从最差的帕累托前沿开始选择待替换个体
    [~, FrontNo] = NSGASort(target_obj);
    max_front = max(FrontNo);

    % 保留 20% 替换比例
    num_replace = ceil(0.2 * pop_size);
    selected_ind = [];

    for front = max_front:-1:1
        inds = find(FrontNo == front);
        selected_ind = [selected_ind; inds(:)];

        if numel(selected_ind) >= num_replace
            selected_ind = selected_ind(1:num_replace);
            break;
        end
    end

    if numel(selected_ind) < num_replace
        remain = setdiff(1:pop_size, selected_ind);
        extra = remain(randperm(numel(remain), ...
            num_replace - numel(selected_ind)));
        selected_ind = [selected_ind; extra(:)];
    end

    % Step 3: 迁移变量，并逐个进行真实评价
    for i = 1:numel(selected_ind)
        ind = selected_ind(i);

        candidate = target_pop(ind, :);
        candidate(selected_idx) = transfer_var;

        exact_obj = fun(candidate);
        assert(isnumeric(exact_obj) && isvector(exact_obj) && ...
               numel(exact_obj) == M, ...
               'STOK:ExactShape', 'fun 必须返回包含 M 个目标值的向量。');

        target_pop(ind, :) = candidate;
        target_obj(ind, :) = reshape(exact_obj, 1, M);
        FEsUsed_add = FEsUsed_add + 1;
    end

    % Step 4: 更新数据库
    database = [target_pop, target_obj];
    success = true;
    source = transfer_source;
end
