
function [child1,child2] = Crossover_Ordered_Operator(parent1, parent2, NDeb)
%CROSSOVER_ORDERED_OPERATOR OX-style crossover for permutation chromosomes.
% Inputs:
%   parent1, parent2 - row/col vectors encoding a permutation (length nvar)
%   NDeb             - max label value (used to repair zeros if any)
% Outputs:
%   child1, child2   - offspring permutations

    nvar = length(parent1);                    % number of genes

    %% ---- child 1 ----
    child1 = zeros(size(parent1));             % preallocate with zeros (unfilled)
    Xpos = sort(randperm(nvar - 1, 2));        % pick 2 cut points within [1..nvar-1]
    child1(Xpos(1):Xpos(2)) = parent1(Xpos(1):Xpos(2));  % copy middle segment from parent1

    j = 1;                                     % pointer into parent2
    for i = 1:nvar
        % fill positions outside the copied segment
        if (i < Xpos(1) || i > Xpos(2)) && j < nvar - 1
            % advance j until parent2(j) not in copied segment
            while sum(parent2(j) == child1(Xpos(1):Xpos(2)))
                j = j + 1;
                if j >= nvar - 1
                    break;
                end
            end
            if j < nvar
                child1(i) = parent2(j);        % place next admissible gene
                j = j + 1;
            end
        end
    end

    % repair any remaining zeros by sampling unused labels
    if any(child1 == 0)
        avOrder = 1:NDeb;                      % candidate labels
        for i = 1:length(child1)
            avOrder(avOrder == child1(i)) = []; % remove labels already used
        end
        child1(child1 == 0) = avOrder(randperm(length(avOrder), sum(child1 == 0)));
    end

    %% ---- child 2 ----
    child2 = zeros(size(parent2));             % preallocate child2
    Xpos = sort(randperm(nvar - 1, 2));        % pick 2 new cut points
    child2(Xpos(1):Xpos(2)) = parent2(Xpos(1):Xpos(2));  % copy middle segment from parent2

    j = 1;                                     % pointer into parent1
    for i = 1:nvar
        if (i < Xpos(1) || i > Xpos(2)) && j < nvar - 1
            
            while sum(parent1(j) == child2(Xpos(1):Xpos(2)))
                j = j + 1;
                if j >= nvar - 1
                    break;
                end
            end
            if j < nvar
                child2(i) = parent1(j);        % place next admissible gene
                j = j + 1;
            end
        end
    end

    % repair zeros if any remain
    if any(child2 == 0)
        avOrder = 1:NDeb;
        for i = 1:length(child2)
            avOrder(avOrder == child2(i)) = [];
        end
        child2(child2 == 0) = avOrder(randperm(length(avOrder), sum(child2 == 0)));
    end
end
