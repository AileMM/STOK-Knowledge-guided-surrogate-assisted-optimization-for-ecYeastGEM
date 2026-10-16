function offspring = Generation(population)
    count = 1;
    MuC =20; MuM =15;
    offspring = zeros(size(population));
 
    N = size(population,1);
    half = fix(N / 2);
    for i = 1:half
        p1 = i;
        p2 = i + half;
        if p2 > N
            p2 = randi(N);  % 或者 p2 = p1; 自交
        end
        % 执行交叉：population(p1,:) 与 population(p2,:)
%       indices = randperm(size(population,1));
%     for i = 1:2:size(population,1)
%         p1 = indices(i);
%         p2 = indices(i+1);
%     for i = 1:ceil(size(population,1) / 2)
%         p1 = i; 
%         p2 = i + fix(size(population,1) / 2);
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

