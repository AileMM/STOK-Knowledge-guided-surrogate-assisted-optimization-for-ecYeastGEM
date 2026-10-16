function [most_similar_idx, mmd_scores] = FindMostSimilarTask(target_pop, source_pops, proj_dim)

% - target_pop: N×D 的目标任务种群（N个个体，D维）
% - source_pops: 一个 cell 数组，每个元素是一个源任务的种群矩阵（N×Di）
% - proj_dim: 投影到的维度（如10）

rng(1);  % 保证随机性可复现
% Step 1: 生成随机投影矩阵 (proj_dim × 原始维度)
D_target = size(target_pop,2);
W = randn(proj_dim, D_target);
% 对目标任务随机投影
target_proj = (W * target_pop')';  % N×proj_dim

num_sources = length(source_pops);
mmd_scores = zeros(num_sources, 1); %保留源任务和目标任务结果

for i = 1:num_sources
    source = source_pops{i};
    D_source = size(source, 2);
 
     % 对源任务进行投影（根据目标任务的维度补零或裁剪）
    if D_source < D_target
        % 源任务维度比目标任务小，补零
        Wi = [source, zeros(size(source,1), D_target - D_source)];
    elseif D_source > D_target
        % 源任务维度比目标任务大，裁剪变量
        Wi = source(:, 1:D_target);
    else
        % 源任务维度等于目标任务维度
        Wi = source;
    end
    source_proj = (W * Wi')';  % N×proj_dim

    % Step 2: 计算 MMD
    mmd_scores(i) = compute_mmd(target_proj, source_proj);
end

[~, most_similar_idx] = min(mmd_scores);
% 检查 most_similar_idx 是否超出范围
if most_similar_idx < 1 || most_similar_idx > num_sources
    error('most_similar_idx 超出了有效范围: %d', most_similar_idx);
end

end

function mmd = compute_mmd(X, Y)
    % 输入： X 和 Y 是 N×d 和 M×d 的矩阵，其中 N 和 M 是样本数量，d 是维度

    % 高斯RBF核函数
    sigma = median([std(X(:)), std(Y(:))]);  % 通过数据计算核宽度
    Kxx = exp(-pdist2(X, X).^2 / (2 * sigma^2));  % N×N 核矩阵
    Kyy = exp(-pdist2(Y, Y).^2 / (2 * sigma^2));  % M×M 核矩阵
    Kxy = exp(-pdist2(X, Y).^2 / (2 * sigma^2));  % N×M 核矩阵

    % 计算 MMD
    mmd = mean(Kxx(:)) + mean(Kyy(:)) - 2 * mean(Kxy(:));
end


