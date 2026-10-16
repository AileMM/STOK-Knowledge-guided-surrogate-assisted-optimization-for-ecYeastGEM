function [best_candidate,best_candidate_obj,paras] =self_opt(database,name_saea,paras,func_surrogate,M,archive_obj)
% 设置不同优化器，框架式改为不同的算法
switch name_saea
    % NSGAII
    case 'NSGAII'
        %% 构建surrogate的时候需要计算相似度开始进行迁移
        
        configs.database = database; %database的实际值
        configs.popsize = 50;
        configs.query_max = 50;
        configs.ini_strategy = 'elite+random';
        configs.acquisition = 'lcb';
        configs.search_engine = 'ea';
        configs.selector = 'roulette_wheel'; 
        
       %% 种群更新--产生子代用代理预测 模型预测值本就是LCB值
        [population,objs_acq] = target_pre(func_surrogate,configs,M);          
        % ================Selection_SPEA2================
        [selected_idx] = Selection_SPEA2(objs_acq, configs.popsize);
        selected_population = population(selected_idx,:);
        objs_acq = objs_acq(selected_idx,:);
        [~,FrontNo] = NSGASort(objs_acq);
        index_offs = FrontNo==1;
        solution_promising = selected_population(index_offs,:);
        objs_promising = objs_acq(index_offs,:);    
        best_candidate = solution_promising;
        best_candidate_obj = objs_promising;

        
end