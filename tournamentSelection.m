function PopulationSelect = tournamentSelection(Population,Ifmax)
%TOURNAMENTSELECTION  Binary tournament selection on scalar fitness.
%   PopulationSelect = TOURNAMENTSELECTION(Population, Ifmax) builds a mating pool
%   the same size as the input by repeated binary tournaments on the .fit field.
%
%   Inputs:
%     Population - struct array with field .fit (scalar fitness)
%     Ifmax      - 1 to keep the higher-fitness candidate, 0 to keep the lower
%   Outputs:
%     PopulationSelect - selected population (same size as input)
PopSize = length(Population);
for n = 1 : PopSize
    s = randi(PopSize,2,1);
    if ~Ifmax

        if Population(s(1)).fit < Population(s(2)).fit
            IDX = s(1);
        else
            IDX = s(2);
        end
    else
        if Population(s(1)).fit < Population(s(2)).fit
            IDX = s(2);
        else
            IDX = s(1);
        end
    end
    PopulationSelect(n) = Population(IDX);
end
end