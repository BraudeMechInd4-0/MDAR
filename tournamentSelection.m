
function PopulationSelect = tournamentSelection(Population, Ifmax)
%TOURNAMENTSELECTION Select individuals using binary tournament.
% Inputs:
%   Population - struct array with field .fit (fitness value)
%   Ifmax      - flag: 1 for maximize fitness, 0 for minimize
% Output:
%   PopulationSelect - selected population (same size as input)

    PopSize = length(Population);          % number of individuals

    %--- Perform tournament for each slot ---
    for n = 1:PopSize
        s = randi(PopSize, 2, 1);          % pick two random candidates

        if ~Ifmax                           % minimizing fitness
            if Population(s(1)).fit < Population(s(2)).fit
                IDX = s(1);
            else
                IDX = s(2);
            end
        else                                % maximizing fitness
            if Population(s(1)).fit < Population(s(2)).fit
                IDX = s(2);
            else
                IDX = s(1);
            end
        end

        PopulationSelect(n) = Population(IDX); % winner goes to next gen
    end
end
