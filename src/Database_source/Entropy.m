function [solution, obj] = Entropy(database, M)

    populations = database(:,1:end-M);
    F1_objs = database(:, end-M+1:end);
    % 归一化目标值到 [0,1]
    F1_objs_norm = (F1_objs - min(F1_objs)) ./ (max(F1_objs) - min(F1_objs));
    % 计算直方图估计概率密度
    numBins = ceil(sqrt(size(F1_objs,1)));  % 经验公式选择分箱数
    [counts, edges] = histcounts(F1_objs_norm, numBins, 'Normalization', 'pdf');
    probs = counts / sum(counts);
    % 计算初始信息熵
    entropy_initial = -sum(probs .* log(probs + eps));
    % 计算每个个体的贡献熵
    contributions = zeros(size(F1_objs,1),1);
    for i = 1:size(F1_objs,1)
        temp_objs = F1_objs_norm;
        temp_objs(i,:) = [];  % 移除个体 i
        % 重新计算概率密度
        [temp_counts, ~] = histcounts(temp_objs, numBins, 'Normalization', 'pdf');
        temp_probs = temp_counts / sum(temp_counts);
        % 计算移除个体后的信息熵
        entropy_removed = -sum(temp_probs .* log(temp_probs + eps));
        % 计算贡献熵（移除后信息熵减少多少）
        % contributions 如果是负的，整体种群的信息熵增加，即该个体对多样性贡献较小，甚至可能是冗余的。
        contributions(i) = entropy_initial - entropy_removed;
    end
%     % 设定熵阈值（保留使熵下降不超过10%的个体）
%     [~, sorted_idx] = sort(abs(contributions), 'descend');
%     selected_idx = sorted_idx(1:topK);
%     threshold = 0.9 * entropy_initial;
%     [~, sorted_idx] = sort(abs(contributions), 'descend');
%     selected_idx = sorted_idx(1:find(cumsum(contributions(sorted_idx)) >= threshold, 1));
    % 选择最终保留下来的个体
    selected_idx = find(contributions > 0);
    solution = populations(selected_idx, :);
    obj = F1_objs(selected_idx, :);

end
