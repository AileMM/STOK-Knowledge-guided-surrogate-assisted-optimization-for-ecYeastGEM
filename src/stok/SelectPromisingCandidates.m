function [best_candidate, best_candidate_obj] = SelectPromisingCandidates(solution_promising, objs_promising, archive_obj, k)
% SelectPromisingCandidates - 从候选解中筛选出非支配且距离当前存档最远的 k 个解
% 输入：
%   solution_promising - 候选解的决策变量 (N × D)
%   objs_promising     - 候选解的预测目标值 (N × M)
%   archive_obj        - 当前存档的真实目标值 (A × M)
%   k                  - 需要选择的候选个体数
% 输出：
%   best_candidate     - 最终选择的 k 个解的决策变量
%   best_candidate_obj - 最终选择的 k 个解的目标值
%   selected_idx       - 选择个体在原始候选集中的索引

    N = size(objs_promising, 1);
    is_dominated = false(N, 1);

    % 判断每个候选解是否被 archive_obj 中的解支配
    for i = 1:N
        for j = 1:size(archive_obj,1)
            if all(archive_obj(j,:) <= objs_promising(i,:)) && any(archive_obj(j,:) < objs_promising(i,:))
                is_dominated(i) = true;
                break;
            end
        end
    end

    % 找出未被支配的个体
    non_dominated_idx = find(~is_dominated);
    best_candidate = solution_promising(non_dominated_idx, :);
    best_candidate_obj = objs_promising(non_dominated_idx, :);

    if isempty(non_dominated_idx)
        % 如果都被支配了，就返回空
        best_candidate = [];
        best_candidate_obj = [];
        selected_idx = [];
        return;
    end
    
    
%     % 计算与 archive 中最小距离作为潜力评分（越大越好）
%     imp = zeros(length(non_dominated_idx), 1);
%     for ii = 1:length(non_dominated_idx)
%         i = non_dominated_idx(ii);
%         dists = vecnorm(archive_obj - objs_promising(i,:), 2, 2);
%         imp(ii) = min(dists);  % 与存档中最近解的距离
%     end
% 
%     % 选择潜力最大的 k 个解
%     [~, sort_idx] = sort(imp, 'descend');
%     top_idx = non_dominated_idx(sort_idx(1:min(k, length(sort_idx))));
% 
%     % 输出对应的解和目标值
%     best_candidate = solution_promising(top_idx, :);
%     best_candidate_obj = objs_promising(top_idx, :);
%     selected_idx = top_idx;
end


