% task_families = {'Sphere','Ellipsoid','Schwefel','Quartic','Ackley','Rastrigin','Griewank','Levy'}; % eight task families
task_families = ...
    {
     'bio1','bio2','bio3','bio4','bio5','bio6','bio7', ...
                            'bio8','bio9','bio10'
    }; 
%这认定目标个数-M
transfer_scenarios = {2,2}; 
% similarity_distributions = {'h1','h2','m1','m2','m3','m4','l1','l2'}; % eight similarity distributions定义了任务之间的相似性分布
k = 10; % the number of source tasks这里对我们代码没有影响应该
specifications = [% 定义源任务：源任务标号、目标个数、维数、源任务数量
    1 2  4261 k; % bio1
    2 2  4261 k; % bio2
    3 1  4261 k; % bio3
    4 2  4261 k; % bio4
    5 1  4261 k; % bio5
    6 2  4261 k; % bio6
    7 1  4261 k; % bio7
    8 2  4261 k; % bio8
    9 2  4261 k; % bio9
    10 2 4261 k;% bio10
    ]; 
folder_problems = '\benchmarks_bio';
sas_names = {'NSGAII','TLRBF','GL-SADE','DDEA-MESS','LSADE','AutoSAEA'}; % six backbone SAS optimizers

num_initial_surrogate = 300; 
num_initial_population = 50; 

FEsMax = 1000; % the maximum number of expensive evaluations
runs = 10; % the number of independent runs

problem_list = 1:size(specifications,1);
num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
no_algorithms = length(sas_names);  %定义不同优化器
% h=waitbar(0,'Starting'); % process monitor
runs_total = no_algorithms*num_problems*runs;
count = 0*no_algorithms*num_problems*runs;

paras.interval = 20; 
paras.ada_flag = 0;

%% Complete Bio1-Bio10 knowledge base, loaded once
        knowledge_base = struct([]);
        for kb_idx = 1:10
            filename = sprintf('bio%d-M2-d4261-k10.mat', kb_idx);  % 文件名
            temp = load(filename, 'knowledge');               % 只加载 knowledge 变量
            data = temp.knowledge;                           % 取出 struct 类型的 knowledge

            if kb_idx == 1
                % 初始化 struct 数组，定义字段结构
                knowledge_base = repmat(data, 1, 10);
            else
                knowledge_base(kb_idx) = data;
            end
        end

% for a =1:no_algorithms
%         stop_tbo = getKnowledgebio(specifications, task_families, transfer_scenarios); 
   for i = 1:length(problem_list)
    % i = 2 ;
    target_idx = problem_list(i);
    fprintf('正在优化第 %d 个任务问题...\n', i);  % <-- 新增这一行
     a = 1;
     results_opt = struct;
     dim = 4261;
     M = 2;
     for r =1:runs 
        % i 循环
%         % 测试函数
%         target_task = problem_family(i,M,dim);
        % bio应用
        target_task = problem_family_bio(target_idx,M,dim);
        name = target_task.name;
        assert(strcmpi(name,sprintf('bio%d',target_idx)), ...
            'STOK:TargetID', 'Target name and ID do not match.');
        fun = target_task.fnc;
        lb = target_task.lb;
        ub = target_task.ub;
        optimum = target_task.F_best;
        rng('shuffle');
        run_rng = rng;   
%         paras.warm_up = lhsdesign(num_initial_surrogate,dim); 
%         paras.warm_up = lb + paras.warm_up.*(ub - lb);    

        %% bio初始化 原模型kcat缩放
        kcats_inf_struct = load('kcats_inf.mat');
        kcats_inf_cell = full(kcats_inf_struct.kcats_inf);  % 转换为普通 cell
        base_kcats = zeros(1, length(kcats_inf_cell));  % 预分配
        for j = 1:length(kcats_inf_cell)
            base_kcats(j) = kcats_inf_cell{j}{3};
        end
        scale_range = [0.001, 10];
        [N, D] = deal(num_initial_population, length(base_kcats));  % M是解个数，D是变量个数
        %%  LHS
        lhs_samples = lhsdesign(N, D);  % 得到 [0,1] 区间的均匀分布样本
        % 将 LHS 样本映射到缩放范围
        scale_factors = scale_range(1) + (scale_range(2) - scale_range(1)) * lhs_samples;
        solutions_initial = scale_factors .* repmat(-base_kcats, N, 1);
%         % ==============================================================================
%         random_scale_factors = scale_range(1) + (scale_range(2) - scale_range(1)) * rand(N, D);
%         solutions_initial = random_scale_factors .* repmat(-base_kcats, N, 1); 
        paras.pop = solutions_initial;
%         paras.pop =  rand(num_initial_population, dim);
%         paras.pop =  lb + paras.pop.*(ub - lb); 
        %% 目标任务的实际优化过程
        % Exclude target only from the source pool; full KB stays intact.
        source_ids = setdiff(1:numel(knowledge_base), target_idx, 'stable');
        candidate_KB = knowledge_base(source_ids);
        assert(numel(source_ids)==9 && ~ismember(target_idx,source_ids), ...
            'STOK:SourcePool', 'Expected nine non-target source tasks.');
        paras.target_id = target_idx;
        paras.source_ids = source_ids;
        target_folder = fullfile('results','Full',sas_names{a});
        if ~exist(target_folder, 'dir')
            mkdir(target_folder);
        end
        paras.fidelity_file = fullfile(target_folder, ...
            sprintf('%s-d%d-FE%d-run%d-fidelity.mat',name,dim,FEsMax,r));
        %% 目标任务的实际优化过程 stop_tbo.knowledge_base
        [database_ktsaea,objs_ktsaea, metric,all_avg_norm,FEs_to_reach,runtime,fidelity]= sas_ckt(fun,lb,ub,FEsMax,sas_names{a},candidate_KB,paras,M,optimum,name);
        results_opt(r).decs = database_ktsaea;
        results_opt(r).objs = objs_ktsaea;
        results_opt(r).igd = metric;
        results_opt(r).all_avg_norm = all_avg_norm;
        results_opt(r).reach = FEs_to_reach;
        results_opt(r).runtime = runtime;
        results_opt(r).variant = 'Full';
        results_opt(r).target_id = target_idx;
        results_opt(r).source_ids = source_ids;
        results_opt(r).run_rng = run_rng;
        results_opt(r).initial_population = solutions_initial;
        results_opt(r).fidelity = fidelity;
        results_opt(r).fidelity_file = paras.fidelity_file;
        results_opt(r).fidelity_metrics = stok_fidelity_metrics(fidelity);
        results_opt(r).FE_final = fidelity.FE_final;
        results_opt(r).FE_record = fidelity.FE_record;
        results_opt(r).norm_obj_record = fidelity.norm_obj_record;
        count = count+1;
        
       % save the results 生成目标文件夹路径
        target_folder = fullfile('results','Full',sas_names{a});
        % 如果目录不存在，则创建
        if ~exist(target_folder, 'dir')
            mkdir(target_folder);
        end
        % 生成目标文件路径
        target_file = fullfile(target_folder, ...
                [task_families{specifications(target_idx,1)},...
                '-d', num2str(specifications(target_idx,3)), '-', num2str(FEsMax), 'k', num2str(specifications(target_idx,4)), '-Fullresults.mat']);
        % 保存文件
        save(target_file, 'results_opt', '-v7.3'); 
     end
  end
% end

function stop_tbo = getKnowledgebio(specifications, task_families, transfer_scenarios)
    persistent cached_knowledge

    if isempty(cached_knowledge)
        cached_knowledge = Knowledge_bio(specifications, task_families, transfer_scenarios);
    end
    stop_tbo = cached_knowledge;
end



% % task_families = {'Sphere','Ellipsoid','Schwefel','Quartic','Ackley','Rastrigin','Griewank','Levy'}; % eight task families
% task_families = ...
%     {
%      'bio1','bio2','bio3','bio4','bio5','bio6','bio7', ...
%                             'bio8','bio9','bio10'
%     }; 
% %这认定目标个数-M
% transfer_scenarios = {2,2}; 
% % similarity_distributions = {'h1','h2','m1','m2','m3','m4','l1','l2'}; % eight similarity distributions定义了任务之间的相似性分布
% k = 10; % the number of source tasks这里对我们代码没有影响应该
% specifications = [% 定义源任务：源任务标号、目标个数、维数、源任务数量
%     1 2  4261 k; % bio1
%     2 2  4261 k; % bio2
%     3 1  4261 k; % bio3
%     4 2  4261 k; % bio4
%     5 1  4261 k; % bio5
%     6 2  4261 k; % bio6
%     7 1  4261 k; % bio7
%     8 2  4261 k; % bio8
%     9 2  4261 k; % bio9
%     10 2 4261 k;% bio10
%     ]; 
% folder_problems = '\benchmarks_bio';
% sas_names = {'NSGAII','TLRBF','GL-SADE','DDEA-MESS','LSADE','AutoSAEA'}; % six backbone SAS optimizers
% 
% num_initial_surrogate = 300; 
% num_initial_population = 50; 
% 
% FEsMax = 1000; % the maximum number of expensive evaluations
% runs = 1; % the number of independent runs
% 
% problem_list = 1:size(specifications,1);
% num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
% no_algorithms = length(sas_names);  %定义不同优化器
% % h=waitbar(0,'Starting'); % process monitor
% runs_total = no_algorithms*num_problems*runs;
% count = 0*no_algorithms*num_problems*runs;
% 
% paras.interval = 20; 
% paras.ada_flag = 0;
% 
% % for a =1:no_algorithms
% %         stop_tbo = getKnowledgebio(specifications, task_families, transfer_scenarios); 
% %    for i =1:length(problem_list)
%     i = 1 ;
%     fprintf('正在优化第 %d 个任务问题...\n', i);  % <-- 新增这一行
%      a = 1;
%      results_opt = struct;
%      dim = 4261;
%      M = 2;
%      for r =1:runs 
%         % i 循环
% %         % 测试函数
% %         target_task = problem_family(i,M,dim);
%         % bio应用
%         target_task = problem_family_bio(i,M,dim);
%         name = target_task.name;
%         fun = target_task.fnc;
%         lb = target_task.lb;
%         ub = target_task.ub;
%         optimum = target_task.F_best;
%         rng('shuffle');   
% %         paras.warm_up = lhsdesign(num_initial_surrogate,dim); 
% %         paras.warm_up = lb + paras.warm_up.*(ub - lb);    
% 
%         %% bio初始化 原模型kcat缩放
%         kcats_inf_struct = load('kcats_inf.mat');
%         kcats_inf_cell = full(kcats_inf_struct.kcats_inf);  % 转换为普通 cell
%         base_kcats = zeros(1, length(kcats_inf_cell));  % 预分配
%         for j = 1:length(kcats_inf_cell)
%             base_kcats(j) = kcats_inf_cell{j}{3};
%         end
%         scale_range = [0.001, 10];
%         [N, D] = deal(num_initial_population, length(base_kcats));  % M是解个数，D是变量个数
%         %%  LHS
%         lhs_samples = lhsdesign(N, D);  % 得到 [0,1] 区间的均匀分布样本
%         % 将 LHS 样本映射到缩放范围
%         scale_factors = scale_range(1) + (scale_range(2) - scale_range(1)) * lhs_samples;
%         solutions_initial = scale_factors .* repmat(-base_kcats, N, 1);
% %         % ==============================================================================
% %         random_scale_factors = scale_range(1) + (scale_range(2) - scale_range(1)) * rand(N, D);
% %         solutions_initial = random_scale_factors .* repmat(-base_kcats, N, 1); 
%         paras.pop = solutions_initial;
% %         paras.pop =  rand(num_initial_population, dim);
% %         paras.pop =  lb + paras.pop.*(ub - lb); 
%         %% 目标任务的实际优化过程
%         knowledge_base = struct([]);
%         for i = 1:10
%             filename = sprintf('bio%d-M2-d4261-k10.mat', i);  % 文件名
%             temp = load(filename, 'knowledge');              % 只加载 knowledge 变量
%             data = temp.knowledge;                           % 取出 struct 类型的 knowledge
% 
%             if i == 1
%                 % 初始化 struct 数组，定义字段结构
%                 knowledge_base = repmat(data, 1, 10);
%             else
%                 knowledge_base(i) = data;
%             end
%         end
%         %% 目标任务的实际优化过程 stop_tbo.knowledge_base
%         [database_ktsaea,objs_ktsaea, metric,all_avg_norm,FEs_to_reach,runtime]= sas_ckt(fun,lb,ub,FEsMax,sas_names{a},knowledge_base,paras,M,optimum,name);
%         results_opt(r).decs = database_ktsaea;
%         results_opt(r).objs = objs_ktsaea;
%         results_opt(r).igd = metric;
%         results_opt(r).all_avg_norm = all_avg_norm;
%         results_opt(r).reach = FEs_to_reach;
%         results_opt(r).runtime = runtime;
%         count = count+1;
% 
%        % save the results 生成目标文件夹路径
%         target_folder = fullfile('.\results\portability\', sas_names{a});
%         % 如果目录不存在，则创建
%         if ~exist(target_folder, 'dir')
%             mkdir(target_folder);
%         end
%         % 生成目标文件路径
%         target_file = fullfile(target_folder, ...
%                 [task_families{specifications(i,1)},...
%                 '-d', num2str(specifications(i,3)), '-1000k', num2str(specifications(i,4)), '-results.mat']);
%         % 保存文件
%         save(target_file, 'results_opt'); 
%      end
% %   end
% % end
% 
% function stop_tbo = getKnowledgebio(specifications, task_families, transfer_scenarios)
%     persistent cached_knowledge
% 
%     if isempty(cached_knowledge)
%         cached_knowledge = Knowledge_bio(specifications, task_families, transfer_scenarios);
%     end
%     stop_tbo = cached_knowledge;
% end
