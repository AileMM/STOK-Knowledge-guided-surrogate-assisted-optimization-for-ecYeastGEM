function source_models = generate_bio_SourceTasks(ecModel_batch, C_sources, N_sources, save_path)
% 生成10个不同碳氮组合的源任务模型，并将其保存为 .mat 文件
% 输入：
%   - ecModel_batch：原始模型
%   - C_sources：碳源 cell 数组
%   - N_sources：氮源 cell 数组
%   - save_path：保存目录（可选，默认当前目录）
% 输出：
%   - source_models：包含10个源任务模型的 cell 数组
% ========修改原模型碳氮源=========
    posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
    posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
    ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),0,'u');
    ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),0,'u');
    if nargin < 4
        save_path = '.';  % 默认保存到当前目录
    end
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    % 固定随机种子确保可复现
    rng(1);
    % 随机选出10个不重复的氮源索引
    N_idx = randperm(length(N_sources), 10);
    % 初始化模型 cell 数组
    source_models = cell(10, 1);

    % 构建并保存每个源任务模型
    for i = 1:10
        C = C_sources{i};          % 第 i 个碳源
        N = N_sources{N_idx(i)};   % 随机选的氮源
        % 设置碳源的摄取通量
        max_uptake = 1000;
        if strcmp(C, '(R)-lactate')
%             vC = [+Inf, +Inf];
            vC = [max_uptake, max_uptake];
        else
%             vC = +Inf;
            vC = max_uptake;
        end
        % 构建模型
        model_i = changeMediaCN(ecModel_batch, C, N, vC, +Inf);
        source_models{i} = model_i;
        % 构造文件名：去除空格和特殊字符
        C_clean = regexprep(C, '[^a-zA-Z0-9]', '_');
        N_clean = regexprep(N, '[^a-zA-Z0-9]', '_');
        filename = fullfile(save_path, ['sourceModel_C_' C_clean '_N_' N_clean '.mat']);
%         % 保存模型变量为 ecModel
%         ecModel = model_i;  
%         save(filename, 'ecModel');
        % 封装为结构体字段
        model_struct = struct();
        model_struct.(['ecModel_C_' C_clean '_N_' N_clean]) = model_i;
        % 保存结构体
        save(filename, '-struct', 'model_struct');
    end
    fprintf('已保存 %d 个源任务模型至目录：%s\n', length(source_models), save_path);
end

function model = changeMediaCN(model,C_source,N_source,C_value,N_value)
    model = openUptake(model,C_source,C_value);
    model = openUptake(model,N_source,N_value);
    if strcmp(C_source,'D-fructose')
        model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
        model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
    elseif strcmp(C_source,'D-mannose')
        model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
        model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
    end
end

function model = openUptake(model,source,value)
    if strcmp(source,'(R)-lactate')
        model = openUptake(model,'(S)-lactate',value(2));
        value = value(1);
    end
    posB = strcmp(model.rxnNames,[source ' exchange (reversible)']);
    posF = strcmp(model.rxnNames,[source ' exchange']);
    % 这个上下界的修改是有问题的--
%     if sum(posB) == 0
%         model = changeRxnBounds(model,model.rxns(posF),-value,'l');
%     else
%         model = changeRxnBounds(model,model.rxns(posB),value,'u');
%     end
%     model = changeRxnBounds(model,model.rxns(posF),0,'u');
    
    if sum(posB) == 0
        % 只有不可逆 exchange
        % 要开启摄取，必须改下界为负值，且上界至少不为0
        model = changeRxnBounds(model, model.rxns(posF), -value, 'l'); % 允许负向摄取
        model = changeRxnBounds(model, model.rxns(posF), value, 'u');   % 保持或打开上界，避免关闭输出
    else
        % 有可逆 exchange
        model = changeRxnBounds(model, model.rxns(posB), -value, 'l'); % 负向摄取允许
        model = changeRxnBounds(model, model.rxns(posB), value, 'u');  % 正向分泌允许
        % 关闭不可逆exchange的上下界，因为可逆exchange已打开摄取和分泌
        model = changeRxnBounds(model, model.rxns(posF), 0, 'l');  % 关闭下界
        model = changeRxnBounds(model, model.rxns(posF), 0, 'u');  % 关闭上界
    end
end


