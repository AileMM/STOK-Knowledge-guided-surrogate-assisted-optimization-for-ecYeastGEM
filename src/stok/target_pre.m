function [population,Obj] = target_pre(func_surrogate,configs,M)

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
    
population = database(:,1:end-M);
Obj = database(:,end-M+1:end);
while no_queries < query_max
    population_offspring = Generation(population);
    objs_acq_offspring = sur_Evaluation(func_surrogate,population_offspring,acquisition,M);  %Ä£ÐÍÔ¤²â
    population = [population; population_offspring];
    Obj = [Obj; objs_acq_offspring];
    [sorted_idx, ~, ~] = NSGASort(Obj);
    population = population(sorted_idx(1:configs.popsize),:);
    Obj = Obj(sorted_idx(1:configs.popsize),:);
    no_queries = no_queries+1;
end


