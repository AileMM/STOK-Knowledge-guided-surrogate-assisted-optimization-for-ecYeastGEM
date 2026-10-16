%% before_after_plot.m
% ---- Step 1: 读取 Excel ----
% 假设迁移前数据在 sheet1，迁移后数据在 sheet2
pre = xlsread('before.xlsx', 1);   % 迁移前 (10×32)
post = xlsread('before.xlsx', 2);  % 迁移后 (10×32)

% ---- Step 2: 清理并展开为向量 ----
preVec  = pre(:);  preVec  = preVec(~isnan(preVec));
postVec = post(:); postVec = postVec(~isnan(postVec));

% ---- Step 3: 箱线图（迁移前） ----
figure;
boxplot(pre, 'Labels', 1:size(pre,2));
xlabel('Var','FontName','Times New Roman','FontSize',16);
ylabel('$k_{cat}$','Interpreter','latex','FontName','Times New Roman','FontSize',16);
title('Before transfer','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',13);
% 简化横坐标（比如只显示每隔4个维度）
grid on; % 打开网格
xticks(1:4:size(pre,2));
print(gcf,'Before_transfer.png','-dpng','-r600');

% ----Step 4: 箱线图（迁移后）----
figure;
boxplot(post, 'Labels', 1:size(post,2));
xlabel('Var','FontName','Times New Roman','FontSize',16);
ylabel('$k_{cat}$','Interpreter','latex','FontName','Times New Roman','FontSize',16);
title('After transfer','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',13);
grid on; % 打开网格
xticks(1:4:size(post,2));
print(gcf,'After_transfer.png','-dpng','-r600');

% ---- Step5: ECDF ----
[f1,x1] = ecdf(preVec);
[f2,x2] = ecdf(postVec);
% KS 检验
[h_ks,p_ks,ksstat] = kstest2(preVec, postVec);

figure; hold on;
h1 = stairs(x1, f1, 'LineWidth',4, 'Color',[0 0.4470 0.7410]);   % 蓝色
h2 = stairs(x2, f2, 'LineWidth',4, 'Color',[0.8500 0.3250 0.0980]); % 橙色
legend([h1,h2], {'Before transfer','After transfer'}, ...
       'Location','best','FontName','Times New Roman','FontSize',13);
xlabel('$k_{cat}$','Interpreter','latex','FontName','Times New Roman','FontSize',16);
ylabel('CDF','FontName','Times New Roman','FontSize',16);
title('CDF of Variables Before and After Transfer', ...
      'FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',13);
grid on; 

% ---- 在图中添加 KS test 结果 ----
xpos = mean(xlim);       % 取横坐标中间
ypos = 0.1;              % 靠下位置
text(xpos, ypos, sprintf('KS test: p=%.2g, stat=%.2f', p_ks, ksstat), ...
     'FontName','Times New Roman','FontSize',14, ...
     'HorizontalAlignment','center');

% ---- Step 6:（可选）保存图像 ----
% 保存为 600 dpi PNG
print(gcf,'ecdf_transfer.png','-dpng','-r600');
