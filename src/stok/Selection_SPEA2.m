function [Next, Fitness] = Selection_SPEA2(Objs, N)
% Environmental selection of SPEA2 without constraints
% Input:
%   - Objs: (M x D) matrix, M individuals with D objective values
%   - N: number of individuals to keep
% Output:
%   - Next: logical index of selected individuals (length M)
%   - Fitness: fitness value for each selected individual

%% Step 1: Calculate fitness
Fitness = CalFitness_no_constraint(Objs);

%% Step 2: Environmental selection
Next = Fitness < 1;
if sum(Next) < N
    [~, Rank] = sort(Fitness);
    Next(Rank(1:N)) = true;
elseif sum(Next) > N
    Del = Truncation(Objs(Next,:), sum(Next) - N);
    Temp = find(Next);
    Next(Temp(Del)) = false;
end

% Output only selected fitness values
Fitness = Fitness(Next);
end

function Del = Truncation(PopObj, K)
% Select part of the solutions by truncation

Distance = pdist2(PopObj, PopObj);
Distance(logical(eye(length(Distance)))) = inf;
Del = false(1, size(PopObj, 1));

while sum(Del) < K
    Remain = find(~Del);
    Temp = sort(Distance(Remain, Remain), 2);
    [~, Rank] = sortrows(Temp);
    Del(Remain(Rank(1))) = true;
end
end


function Fitness = CalFitness_no_constraint(PopObj)
% Calculate fitness for unconstrained multi-objective optimization
N = size(PopObj, 1);
%% Dominance Matrix
Dominate = false(N);
for i = 1:N-1
    for j = i+1:N
        k = any(PopObj(i,:) < PopObj(j,:)) - any(PopObj(i,:) > PopObj(j,:));
        if k == 1
            Dominate(i,j) = true;
        elseif k == -1
            Dominate(j,i) = true;
        end
    end
end

%% Strength S(i)
S = sum(Dominate, 2);

%% Raw fitness R(i)
R = zeros(1, N);
for i = 1:N
    R(i) = sum(S(Dominate(:, i)));
end

%% Density D(i)
Distance = pdist2(PopObj, PopObj);
Distance(logical(eye(N))) = inf;
Distance = sort(Distance, 2);
D = 1 ./ (Distance(:, floor(sqrt(N))) + 2);

%% Final Fitness
Fitness = R + D';
end

