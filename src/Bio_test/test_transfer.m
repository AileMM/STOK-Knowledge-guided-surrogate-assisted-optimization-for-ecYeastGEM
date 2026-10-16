%% 读取数据
% 假设文件名为 data.xlsx
% 两行数据：第1行 -> 向量1, 第2行 -> 向量2
data = readmatrix('data.xlsx');

x1 = data(1, :)';  % 转置成列向量
x2 = data(2, :)';

%% 1. 绘制累积分布函数 (ECDF)
figure;
hold on;
ecdf(x1);
ecdf(x2);
legend('Before transfer','After transfer','Location','best');
xlabel('%k_{cat}%');
ylabel('CDF');
title('Compare the %k_{cat}% with knowledge transfer');
grid on;

%% 2. Kolmogorov-Smirnov 检验
[h,p,ksstat] = kstest2(x1, x2);
fprintf('KS 检验结果: h=%d, p=%.6f, KS统计量=%.4f\n', h, p, ksstat);

% h = 1 表示两个向量分布有显著差异
% h = 0 表示不能拒绝相同分布的假设

%% 3. 额外差异度量
% 均值差异
mean_diff = mean(x1) - mean(x2);

% KL 散度 (基于直方图离散化)
edges = linspace(min([x1;x2]), max([x1;x2]), 100); % 100个区间
p1 = histcounts(x1, edges, 'Normalization', 'probability');
p2 = histcounts(x2, edges, 'Normalization', 'probability');
kl_div = sum(p1 .* log((p1+eps)./(p2+eps))); % 避免除零

fprintf('均值差: %.4f, KL散度: %.4f\n', mean_diff, kl_div);

%% 4. 可视化直方图比较
figure;
histogram(x1, edges, 'Normalization','probability', 'FaceAlpha',0.5);
hold on;
histogram(x2, edges, 'Normalization','probability', 'FaceAlpha',0.5);
legend('Before transfer','After transfer','Location','best');
xlabel('%k_{cat}%');
ylabel('Probability');
title('Histogram of the distributions of two vectors');
grid on;
