%% 用于算法的迭代曲线，把十个目标任务结果合并到一起
% 设置目标文件夹路径 
targetFolder = 'C:\Users\ADcalcium\Desktop\Prior knowledge\代码\MTO-new\results\portability\NSGAII';

% 只获取符合命名规则的 mat 文件（如 bio1-xxxx-results.mat）
matFiles = dir(fullfile(targetFolder, 'bio*-results.mat'));
numFiles = length(matFiles);

allData = [];              % 初始化数据矩阵
fileNames = {};            % 用于存储文件名

for i = 1:numFiles
    filePath = fullfile(matFiles(i).folder, matFiles(i).name);
    data = load(filePath);  % 加载 mat 文件
    
    % 获取结构体名称（假设只有一个主结构体）
    fieldnamesData = fieldnames(data);
    mainStruct = data.(fieldnamesData{1});

    % 提取目标字段
    if isfield(data, 'results_opt') && isfield(data.results_opt, 'all_avg_norm')
        vec = data.results_opt.all_avg_norm;
        allData(i, 1:length(vec)) = vec;  % 自动补 NaN（向右扩展）
        fileNames{i,1} = matFiles(i).name;  % 存储文件名
    else
        error(['Missing field in: ', matFiles(i).name]);
    end
end

% 构建表格并加上文件名列
T = array2table(allData);
T = addvars(T, fileNames, 'Before', 1, 'NewVariableNames', 'FileName');

% 输出到 Excel
outputFile = 'combined_results.xlsx';
writetable(T, outputFile);

disp(['? 成功写入 ', outputFile]);
