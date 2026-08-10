function [Population, BestSolAfter] = EliteProcedure(OffspringPop,Population,BestSol,GenParam)
%ELITEPROCEDURE  Elitist environmental selection for the single-objective GA.
%   [Population,BestSolAfter] = ELITEPROCEDURE(OffspringPop, Population, BestSol, GenParam)
%   combines parents, offspring, and the incumbent best, then keeps GenParam.PopSize
%   individuals for the next generation according to GenParam.EliteMode.
%
%   Inputs:
%     OffspringPop - struct array of offspring individuals
%     Population   - current-generation population
%     BestSol      - best individual from the previous generation
%     GenParam     - struct with fields:
%                      .PopSize   - target population size
%                      .EliteMode - "all" | "none" | "One" | other (default "H")
%   Outputs:
%     Population    - next-generation population after elitism
%     BestSolAfter  - best individual across the combined pool (for logging)
%
%   EliteMode: "all" selects from [Population, BestSol, Offspring]; "none" uses
%   Offspring as-is; "One" lets BestSol compete with Offspring only; any other value
%   keeps the front half of the current population plus BestSol and Offspring.
PopSize = GenParam.PopSize;
EliteMode = GenParam.EliteMode;
PopAll = [Population,BestSol,OffspringPop];
Fall = [PopAll(:).fit];
[~,idx] = sort(Fall);
BestSolAfter = PopAll(idx(1));
%save('temp','BestSolAfter');
if strcmp(EliteMode, "all")
    Population = PopAll(idx(1:PopSize));
elseif strcmp(EliteMode, "none")
    Population = OffspringPop;
elseif strcmp(EliteMode, "One")
    PopAll = [BestSol,OffspringPop];
    Fall = [PopAll(:).fit];
    [~,idx] = sort(Fall);
    Population = PopAll(idx(1:PopSize));
else
    PopAll = [Population(2:floor(length(Population)/2)),BestSol,OffspringPop];
    Fall = [PopAll(:).fit];
    [~,idx] = sort(Fall);
    Population = PopAll(idx(1:PopSize));
end
end