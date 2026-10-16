% function [sorted_idx, similarity_scores] = compare_PFs_by_cosine(target_PF, source_PF_list)
% % target_PF: N_target × M 矩阵，目标任务的预测 PF
% % source_PF_list: 1 × K cell，每个 cell 是一个 N_source × M 的源任务 PF
% % sorted_idx: 相似度从高到低排序的索引
% % similarity_scores: 每个源任务的余弦相似度分数
% 
%     % 计算目标 PF 的中心点向量
%     center_target = mean(target_PF, 1);
%     % 归一化目标中心向量
%     center_target = center_target / norm(center_target);
%     
%     K = length(source_PF_list);
%     similarity_scores = zeros(K,1);
%     for i = 1:K
%         source_PF = source_PF_list{i};
%         center_source = mean(source_PF.objs, 1);
%         center_source = center_source / norm(center_source);     
%         % 计算余弦相似度（向量点积）
%         similarity_scores(i) = dot(center_target, center_source);
%     end
%     
%     % 按照余弦相似度降序排序
%     [~, idx] = sort(similarity_scores, 'descend');
%     sorted_idx = idx(1);
%     sorted_idx = idx;
% end

function [sorted_idx, similarity_scores] = compare_PFs_by_cosine(target_PF, source_PF_list)
    target_PF = normalize(target_PF, 2, 'norm');  % 每一行归一化
    K = length(source_PF_list);
    similarity_scores = zeros(K,1);
    
    for i = 1:K
        source_PF = source_PF_list{i}.objs;
        source_PF = normalize(source_PF, 2, 'norm');  % 每一行归一化
        
        % 所有点对之间的余弦相似度矩阵
        sim_matrix = target_PF * source_PF';
        max_sim = max(sim_matrix, [], 2);  % 每个目标点与源PF中最相似的点
        similarity_scores(i) = mean(max_sim);  % 平均最大相似度
    end
    
    [~, idx] = sort(similarity_scores, 'descend');
    sorted_idx = idx(1);
end

