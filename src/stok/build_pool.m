% task_families = {'Sphere','Ellipsoid','Schwefel','Quartic','Ackley','Rastrigin','Griewank','Levy'}; % eight task families
task_families = ...
    {
    'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7',...
    'WFG1','WFG2','WFG3','WFG4',...
    'WFG5',...
    'ZDT1','ZDT2','ZDT3',...
    'ZDT4','ZDT6'
    }; 
%这认定目标个数-M
transfer_scenarios = {2,3}; 
% similarity_distributions = {'h1','h2','m1','m2','m3','m4','l1','l2'}; % eight similarity distributions定义了任务之间的相似性分布
k = 17; % the number of source tasks这里对我们代码没有影响应该
specifications = [% 定义源任务：源任务标号、目标个数、维数、源任务数量
    1 2  50 k; % DTLZ1
    2 2  25 k; % DTLZ2
    3 1  30 k; % DTLZ3
    4 2  30 k; % DTLZ4
    5 1  25 k; % DTLZ5
    6 2  30 k; % DTLZ6
    7 1  25 k; % DTLZ7
    8 2  30 k; % WFG1
    9 2  25 k; % WFG2
    10 2  12 k; % WFG3
    11 2  50 k; % WFG4
    12 1  30 k; % WFG5
%     13 1  35 k; % ZDT1
%     14 1  30 k; % ZDT2
%     15 1  50 k; % ZDT3
%     16 1  60 k; % ZDT4
%     17 1  50 k; % ZDT6
    ]; 
folder_problems = '\benchmarks';
sas_names = {'NSGAII','TLRBF','GL-SADE','DDEA-MESS','LSADE','AutoSAEA'}; % six backbone SAS optimizers

num_initial_surrogate = 300; 
num_initial_population = 50; 

FEsMax =1000; % the maximum number of expensive evaluations
runs = 3; % the number of independent runs

problem_list = 1:size(specifications,1);
num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
no_algorithms = length(sas_names);  %定义不同优化器
% h=waitbar(0,'Starting'); % process monitor
runs_total = no_algorithms*num_problems*runs;
count = 0*no_algorithms*num_problems*runs;

paras.interval = 20; 
paras.ada_flag = 0;
%% 重构知识库
% for a =1:no_algorithms
%         stop_tbo = getKnowledge(specifications, task_families, transfer_scenarios); 
%% 目标任务求解
   for i =1:length(problem_list)
%    for i = 11
     fprintf('正在优化第 %d 个任务问题...\n', i);  % <-- 新增这一行
     a = 1;
     results_opt = struct;
     dim = 10;
     M = 2;
     
    all_decs = [];
    all_objs = [];
    all_igd  = [];
    run_ids  = [];
    
     for r =1:runs 
        % i 循环
        target_task = problem_family(i,M,dim);
        fun = target_task.fnc;
        lb = target_task.lb;
        ub = target_task.ub;
        optimum = target_task.F_best;
        rng('shuffle');   
        paras.warm_up = lhsdesign(num_initial_surrogate,dim); 
        paras.warm_up = lb + paras.warm_up.*(ub - lb);    
        paras.pop =  rand(num_initial_population, dim);
        paras.pop =  lb + paras.pop.*(ub - lb);
        
        %% 目标任务的实际优化过程
        [database_ktsaea,objs_ktsaea, metric]= sas_ckt(fun,lb,ub,FEsMax,sas_names{a},stop_tbo.knowledge_base,paras,M,optimum);
        results_opt(r).decs = database_ktsaea;
        results_opt(r).objs = objs_ktsaea;
        results_opt(r).igd = metric;
        
%         n_i = size(database_ktsaea, 1);
%         all_decs = [all_decs; database_ktsaea];
%         all_objs = [all_objs; objs_ktsaea];
%         all_igd  = [all_igd; metric];
%         run_ids  = [run_ids; r*ones(n_i,1)];  % 每一行标记所属 run
     end
         
       % save the results 生成目标文件夹路径
        target_folder = fullfile('.\results\portability\', sas_names{a});
        % 如果目录不存在，则创建
        if ~exist(target_folder, 'dir')
            mkdir(target_folder);
        end
        % 生成目标文件路径
        target_file = fullfile(target_folder, ...
                [task_families{specifications(i,1)}, '-M', num2str(specifications(i,2)), '-no-transfer',...
                '-d', num2str(specifications(i,3)), '-1000k', num2str(specifications(i,4)), '-results.mat']);
        % 保存文件
        save(target_file, 'results_opt'); 

%      end
  end
% end
function stop_tbo = getKnowledge(specifications, task_families, transfer_scenarios)
    persistent cached_knowledge

    if isempty(cached_knowledge)
        cached_knowledge = Knowledge(specifications, task_families, transfer_scenarios);
    end
    stop_tbo = cached_knowledge;
end
