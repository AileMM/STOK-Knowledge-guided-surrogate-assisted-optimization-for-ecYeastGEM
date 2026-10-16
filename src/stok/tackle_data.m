%% 读取已经构建好的知识库加入  设定包含 mat 文件的文件夹路径
% 设置 .mat 文件目录
folderPath = 'C:\Users\ADcalcium\Desktop\Prior knowledge\代码\MTO-new\benchmarks_bio';
matFiles = dir(fullfile(folderPath, '*.mat'));

knowledge_base = [];  % 初始化为空 struct 数组

for i = 1:length(matFiles)
    filePath = fullfile(matFiles(i).folder, matFiles(i).name);
    data = load(filePath);  % 加载当前 .mat 文件
    
    if isfield(data, 'knowledge')
        if isempty(knowledge_base)
            knowledge_base = data.knowledge;  % 第一次直接赋值
        else
            knowledge_base(end+1) = data.knowledge;  % 后续结构数组追加
        end
    else
        warning(['Missing "knowledge" field in ', matFiles(i).name]);
    end
end
% 存入 base 工作区，供后续调用
assignin('base', 'knowledge_base', knowledge_base);
disp(['? 成功加载 ', num2str(length(knowledge_base)), ' 个 knowledge struct']);

% 
% %% ========================================================================
% inputFolder = 'C:\Users\ADcalcium\Desktop\Prior knowledge\代码\MTO-new\results\portability\NSGAII';     % <== 鏇挎崲涓轰綘鐨? .mat 鏂囦欢澶?
% outputFolder = 'C:\Users\ADcalcium\Desktop\Prior knowledge\代码\MTO-new\results\portability\NSGAII';  % <== 鏇挎崲涓鸿緭鍑? .xlsx 鏂囦欢澶?
% 
% if ~exist(outputFolder, 'dir')
%     mkdir(outputFolder);
% end
% 
% % 鑾峰彇鎵?鏈? .mat 鏂囦欢
% matFiles = dir(fullfile(inputFolder, 'bio1*.mat'));
% 
% for i = 1:length(matFiles)
% 
%     matPath = fullfile(inputFolder, matFiles(i).name);
%     [~, baseFileName, ~] = fileparts(matFiles(i).name);
%     excelPath = fullfile(outputFolder, [baseFileName, '.xlsx']);
%     
%     % 鍔犺浇 struct 鏁版嵁
%     S = load(matPath);
%     fieldName = fieldnames(S);
%     data = S.(fieldName{1}); 
%     
%     writematrix(data.decs,         excelPath, 'Sheet', 'decs');
%     writematrix(data.objs,         excelPath, 'Sheet', 'objs');
%     writematrix(data.igd(:),       excelPath, 'Sheet', 'igd'); 
%     writematrix(data.all_avg_norm, excelPath, 'Sheet', 'all_avg_norm');
%     
%     % 提取 reach 和 runtime
%     reach = data.reach;
%     runtime = data.runtime;
% 
%     % 如果为空，用 NaN 占位
%     if isempty(reach)
%         reach = NaN;
%     end
%     if isempty(runtime)
%         runtime = NaN;
%     end
% 
%     % 构建 meta_info，其中第二行为 reach 和 runtime 数据
%     meta_info = {
%         "reach",   "runtime";
%         reach,     runtime
%     };
%     writecell(meta_info, excelPath, 'Sheet', 'meta_info');
%     
%     fprintf('宸插鐞? %s 骞跺鍑轰负 Excel銆俓n', baseFileName);
% end

