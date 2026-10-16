classdef Knowledge
    properties
        % 配置参数
        trans_sce;  % 迁移场景
        dim;  % 任务维度
        k;  % 源任务个数
        optimizer = 'saea_source';  % 优化算法
        FEsMax = 100;  % 最大评估次数
        mode = 'opt';  % 优化模式
        folder_stops = 'benchmarks'; % 知识库文件夹
        knowledge_base = struct;  % 知识库
        kb_idx
        total_tasks
    end
    properties(SetAccess = protected)
        % 任务集合
        problem_families = {
                            'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7', ...
                            'WFG1','WFG2','WFG3','WFG4',...
                            'WFG5', ...
                            'ZDT1','ZDT2','ZDT3',...
                            'ZDT4','ZDT6'
                            };
    end
    methods
        function obj = Knowledge(specifications, task_families, transfer_scenarios)
            for t = 1:size(specifications, 1)
                spec = specifications(t, :);  % 获取第t个任务的配置
                obj.trans_sce = transfer_scenarios{spec(2)};  % 迁移场景
                obj.dim = spec(3);  % 任务维度
                obj.k = spec(4);  % 每个任务的优化次数
                % 调用 Generation 生成源任务并优化
                obj = obj.Generation(task_families, spec, t);  % 传递任务配置 spec 和任务索引 t
            end
        end
        function obj = Generation(obj, task_families, spec, t)
           source_task = task_families{spec(1)};  
            % 调用 problem_family 函数生成问题
            problem = problem_family(find(strcmp(obj.problem_families, source_task)), obj.trans_sce, obj.dim);
            for i = 1:obj.k
                % 优化过程
                cmd_opt = ['[solutions, objs, surrogate, selected_variables, selected_C, selected_D] = ', obj.optimizer, '(problem, obj.FEsMax, obj.trans_sce, i);'];
                eval(cmd_opt);  % saea_source

                % 存储优化结果
                obj.knowledge_base(t).solutions = solutions;
                obj.knowledge_base(t).objs = objs;
                obj.knowledge_base(t).surr = surrogate;
                obj.knowledge_base(t).key_var = selected_variables;
                obj.knowledge_base(t).var_c = selected_C;
                obj.knowledge_base(t).var_d = selected_D;

                fprintf('BO-LCB is solving the source task: %s\n', source_task);
            end
            knowledge = obj.knowledge_base(t);
            % 保存结果到文件
            save_filename = fullfile(obj.folder_stops, ...
                [source_task, '-M', num2str(obj.trans_sce),  ...
                 '-d', num2str(obj.dim), '-k', num2str(obj.k), '.mat']);
             
            save(save_filename, 'solutions', 'objs', 'knowledge');
        end

    end
end

