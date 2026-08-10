function MatingPool = SelectionByRank(Population)
%SELECTIONBYRANK  NSGA-II binary tournament selection (rank, then crowding).
%   MatingPool = SELECTIONBYRANK(Population) builds a mating pool the same size as
%   the input by binary tournaments: the lower non-domination rank wins; ties are
%   broken by the larger crowding distance.
%
%   Inputs:
%     Population - struct array with fields .Rank and .CrowdDis
%   Outputs:
%     MatingPool - selected population (same size as input)
N = length(Population);
MatingPool(N) = Population(1);
idx = randi(N,N,2);
for i = 1 : N
    if Population(idx(i,1)).Rank < Population(idx(i,2)).Rank
        MatingPool(i) = Population(idx(i,1));
    elseif Population(idx(i,1)).Rank > Population(idx(i,2)).Rank
        MatingPool(i) = Population(idx(i,2));
    elseif Population(idx(i,1)).CrowdDis > Population(idx(i,2)).CrowdDis
        MatingPool(i) = Population(idx(i,1));
    else
        MatingPool(i) = Population(idx(i,2));
    end
end
