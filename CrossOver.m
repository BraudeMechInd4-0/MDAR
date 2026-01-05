
function PopulationAfterCross = CrossOver(Population, GenParam)
%CROSSOVER Apply crossover on pairs: SBX for waits + OX for order.
% Inputs:
%   Population   - struct array with fields .WaitUntilManuver, .Order
%   GenParam     - struct with fields:
%                  .pc    -> crossover probability (0..1)
%                  .NDeb  -> number of debris labels
%                  .PopSize -> population size
%                  .T1max -> upper bound for wait times (lower bound = 1 here)
% Output:
%   PopulationAfterCross - population after crossover

    pc   = GenParam.pc;                 % crossover probability
    NDeb = GenParam.NDeb;               % label count for permutation repair
    PopSize = GenParam.PopSize;         % population size
    T1limits = [1 GenParam.T1max];      % SBX bounds for wait times

    PopulationAfterCross = Population;  % initialize output

    % Pairwise crossover: (1,2), (3,4), ...
    for i = 1:2:PopSize
        % With prob > pc perform crossover; else copy parents (as given)
        if rand() > pc
            % Simulated Binary Crossover on wait vectors (real-valued)
            [PopulationAfterCross(i).WaitUntilManuver, ...
             PopulationAfterCross(i+1).WaitUntilManuver] = ...
                SBX(Population(i).WaitUntilManuver, ...
                    Population(i+1).WaitUntilManuver, 1, T1limits);

            % Ordered crossover on permutation chromosomes (integer labels)
            [PopulationAfterCross(i).Order, ...
             PopulationAfterCross(i+1).Order] = ...
                Crossover_Ordered_Operator(Population(i).Order, ...
                                           Population(i+1).Order, NDeb);
        else
            % No crossover: children are clones of parents
            PopulationAfterCross(i).WaitUntilManuver   = Population(i).WaitUntilManuver;
            PopulationAfterCross(i+1).WaitUntilManuver = Population(i+1).WaitUntilManuver;
            PopulationAfterCross(i).Order              = Population(i).Order;
            PopulationAfterCross(i+1).Order            = Population(i+1).Order;
        end
    end
end
