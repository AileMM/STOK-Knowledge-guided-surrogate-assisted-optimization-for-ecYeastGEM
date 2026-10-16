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
    13 1  35 k; % ZDT1
    14 1  30 k; % ZDT2
    15 1  50 k; % ZDT3
    16 1  60 k; % ZDT4
    17 1  50 k; % ZDT6
    ]; 
folder_problems = '.\benchmarks';
sas_names = {'BO-LCB','TLRBF','GL-SADE','DDEA-MESS','LSADE','AutoSAEA'}; % six backbone SAS optimizers
num_initial_solutions = 50; % initial solutions POP
FEsMax = 100; % the maximum number of expensive evaluations
problem_list = 1:size(specifications,1);
num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
no_algorithms = length(sas_names);  %定义不同优化器
% h=waitbar(0,'Starting'); % process monitor
runs = 2; % the number of independent runs

runs_total = no_algorithms*num_problems*runs;
count = 0*no_algorithms*num_problems*runs;
paras.interval = 20; % the transfer interval in terms of the true evaluation
paras.ada_flag = 0; % whether to use the task adaptation, 1: yes, 0: no

runs = 2; % the number of independent runs
for a =1:no_algorithms
    for i =1:length(problem_list)
        d = [10,20,30,50,100];
        dim = d(randi(5));
        num = [2, 3];
        M = num(randi(2));
        target_task = problem_family(i,M,dim);
        fun = target_task.fnc;
        lb = target_task.lb;
        ub = target_task.ub;
        results_opt = struct;
        paras.warm_up = lhsdesign_modified(num_initial_solutions,lb,ub); % 目标任务的初始化种群
         [database_ktsaea,objs_ktsaea,states_ktsaea]= sas_ckt(fun,lb,ub,FEsMax,...
        sas_names{a},stop_tbo.knowledge_base,paras,M);
        results_opt(r).objs = objs_ktsaea;
        results_opt(r).states = states_ktsaea;
        count = count+1;

    end
%         for r = 1:runs  
%            
%         end
%         
        % save the results 生成目标文件夹路径
        target_folder = fullfile('.\results\portability\', sas_names{a});

        % 如果目录不存在，则创建
        if ~exist(target_folder, 'dir')
            mkdir(target_folder);
        end

        % 生成目标文件路径
        target_file = fullfile(target_folder, ...
            [task_families{specifications(n,1)}, '-T', transfer_scenarios{specifications(n,2)}, ...
            '-S', similarity_distributions{specifications(n,3)}, '-d', num2str(specifications(n,4)), ...
            '-k', num2str(specifications(n,5)), '-results.mat']);

        % 保存文件
        save(target_file, 'results_opt');
        
    
end