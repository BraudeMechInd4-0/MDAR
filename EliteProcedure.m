
function [Population, BestSolAfter] = EliteProcedure(OffspringPop, Population, BestSol, GenParam)
%ELITEPROCEDURE Combine parents/offspring and keep elites for next gen.
% Inputs:
%   OffspringPop - struct array of offspring individuals
%   Population   - current generation population
%   BestSol      - best individual from previous generation
%   GenParam     - struct with fields:
%                  .PopSize   -> target population size
%                  .EliteMode -> 'all' | 'none' | 'One' | 'H' (fallback/other)
% Outputs:
%   Population   - next generation population after elitism
%   BestSolAfter - best individual among combined pools (for logging)
%
% Notes on EliteMode:
%   'all' : select best PopSize from [Population, BestSol, OffspringPop]
%   'none': ignore elitism; use OffspringPop as-is
%   'One' : ensure BestSol competes with Offspring only
%   other : keep half of current Population (excluding first), + BestSol + Offspring,
%           then select best PopSize (custom "H" behavior)

    PopSize  = GenParam.PopSize;               % Target population size
    EliteMode = GenParam.EliteMode;            % Elitism strategy

    %--- Combine and rank by fitness (descending = higher is better) ---
    PopAll = [Population, BestSol, OffspringPop];   % Pool all candidates
    Fall   = [PopAll(:).fit];                        % Collect fitness values
    [~, idx] = sort(Fall, 'descend');                % Rank by fitness
    BestSolAfter = PopAll(idx(1));                   % Keep global best

    %--- Apply elitism strategy to form next population ---
    if strcmp(EliteMode, "all")
        % Keep top PopSize from full pool
        Population = PopAll(idx(1:PopSize));

    elseif strcmp(EliteMode, "none")
        % No elitism; next gen is only offspring
        Population = OffspringPop;

    elseif strcmp(EliteMode, "One")
        % BestSol competes only with offspring; take top PopSize
        PopAll = [BestSol, OffspringPop];
        Fall   = [PopAll(:).fit];
        [~, idx] = sort(Fall, 'descend');
        Population = PopAll(idx(1:PopSize));

    else
        % Custom "H" style: keep mid slice of current population + BestSol + Offspring
        % (drops the first member, takes up to ~half for diversity)
        midStop = floor(length(Population)/2);
        baseSeg = Population(2:max(2, midStop));     % safe slice if Pop is small
        PopAll  = [baseSeg, BestSol, OffspringPop];
        Fall    = [PopAll(:).fit];
        [~, idx] = sort(Fall, 'descend');
        Population = PopAll(idx(1:PopSize));
    end
end
