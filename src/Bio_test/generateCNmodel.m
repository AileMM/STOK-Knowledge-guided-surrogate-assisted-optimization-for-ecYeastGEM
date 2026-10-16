%% 生成不同碳氮源模型
load('ecYeastGEM_batch.mat');   % ecModel_batch
load('C_sources.mat');          % C_sources
load('N_sources.mat');          % N_sources
% 指定保存路径，例如 'source_models'
% mkdir('bio_source_models');
source_models = generate_bio_SourceTasks(ecModel_batch, C_sources, N_sources, 'bio_source_models');

% %% 检查模型碳氮源是否修改成功
% data = load('sourceModel_C__R__lactate_N_L_aspartate.mat');
% vars = fieldnames(data);       % 获取结构体中的所有变量名
% ecModel = data.(vars{1});      % 取第一个变量的值赋给 ecModel
% checkMediaChange(ecModel, '(R)-lactate', 'L-aspartate');
% 
% % sol = optimizeCbModel(ecModel, 'max');