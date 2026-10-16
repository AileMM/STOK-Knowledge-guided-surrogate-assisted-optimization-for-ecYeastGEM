function [new_population,Obj] = transfer_pre(func_surrogate,configs,M,transfer_source)
    % initialization
    database = configs.database;
    popsize = configs.popsize;
    query_max = configs.query_max;
    ini_strategy = configs.ini_strategy;
    acquisition = configs.acquisition;
    search_engine = configs.search_engine;
    selector = configs.selector;
    FEsUsed = size(database,1);
    no_queries = 0;
    % 目标任务x和y
    population = database(:,1:end-M);
    Obj = database(:,end-M+1:end);
    
    x_transfer = transfer_source.solutions;
    D_target = size(population, 2);
    D_source = size(x_transfer, 2);
    if D_source > D_target
        X_transferred = x_transfer(:, 1:D_target);
    elseif D_source < D_target
        % 随机补全剩余维度
        n_ind = size(x_transfer, 1);
        X_transferred = [x_transfer, rand(n_ind, D_target - D_source)];
    else
        X_transferred = x_transfer;
    end
    
   new_population = [population;X_transferred];
    while no_queries < query_max
        population_offspring = Generation(new_population);
        objs_acq_offspring = sur_Evaluation(func_surrogate,population_offspring,acquisition,M);  
        new_population = [new_population; population_offspring];
        Obj = [Obj; objs_acq_offspring]; 
        % ======这里修改一下===========
        [rank,~] = NSGASort(Obj);
        new_population = new_population(rank(1:configs.popsize),:);
        Obj = Obj(rank(1:configs.popsize),:);
        no_queries = no_queries+1;
    end
end

function offspring = Generation(population)
    count = 1;
    MuC =20; MuM =15;
    offspring = zeros(size(population));
    for i = 1:ceil(size(population,1) / 2)
        p1 = i; 
        p2 = i + fix(size(population,1) / 2);
%         p2 = mod(i + fix(length(population) / 2) - 1, length(population)) + 1;
        offspring(count,:) = population(p1,:);
        offspring(count + 1,:) = population(p2,:);
        [offspring(count,:), offspring(count + 1,:)] = GA_Crossover(population(p1,:), population(p2,:), MuC);
        offspring(count,:) = GA_Mutation(offspring(count,:), MuM);
        offspring(count + 1,:) = GA_Mutation(offspring(count + 1,:), MuM);
        offspring(count,:) = min(max(offspring(count,:), 0), 1);
        offspring(count + 1,:) = min(max(offspring(count + 1,:), 0), 1);
        count = count + 2;  % 更新索引
    end
end





