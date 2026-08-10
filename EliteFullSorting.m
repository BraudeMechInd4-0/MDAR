function PopulationUpdated = EliteFullSorting(Population,Offsprings,N,elitePer)
%ELITEFULLSORTING  NSGA-II environmental selection (rank + crowding truncation).
%   PopulationUpdated = ELITEFULLSORTING(Population, Offsprings, N, elitePer) merges
%   the current non-dominated set with the offspring, recomputes non-domination rank
%   and crowding distance (CalcRankAndDistance), and keeps the best N individuals
%   ordered by ascending rank then descending crowding distance.
%
%   Inputs:
%     Population - current elite / non-dominated struct array
%     Offsprings - offspring struct array to merge in
%     N          - target population size
%     elitePer   - elitism fraction; 1 keeps the best N by rank/crowding, otherwise
%                  fills (1-elitePer)*N from offspring and the remainder at random
%   Outputs:
%     PopulationUpdated - next-generation population of size N
Population = CalcRankAndDistance(Population);
AllPop = [Population,Offsprings];
AllPop = CalcRankAndDistance(AllPop);
AllPop(isnan([AllPop(:).CrowdDis]))=[];
RankAndDistMatrix = [[AllPop(:).Rank]',[AllPop(:).CrowdDis]'];
[~, I ] = sortrows(RankAndDistMatrix,[1 -2]);

if elitePer == 1

    if length(I)>N
        PopulationUpdated = AllPop(I(1:N));
    else
        PopulationUpdated = AllPop;
        while length(PopulationUpdated) < N
            PopulationUpdated = [PopulationUpdated,AllPop(randi(length(AllPop)))];
        end
    end
else


    PopulationUpdated = Offsprings(1:round((1-elitePer)*N));

    while length(PopulationUpdated) < N
        PopulationUpdated = [PopulationUpdated,AllPop(randi(length(AllPop)))];
    end
end

