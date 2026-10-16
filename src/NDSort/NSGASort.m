function [rank, FrontNo, CrowdDis] = NSGASort(population)
FrontNo = NDSort(population,inf);
CrowdDis = CrowdingDistance(population, FrontNo);
[~, rank] = sortrows([FrontNo', -CrowdDis']);
end
