
function [offspring1, offspring2] = SBX(parent1, parent2, CrossOverRate, Limits)
%SBX Simulated Binary Crossover for real-valued vectors.
% Inputs:
%   parent1, parent2  - real-valued row/column vectors (same length)
%   CrossOverRate     - probability of applying SBX (0..1)
%   Limits            - [yl yu] lower/upper bounds for each gene
% Outputs:
%   offspring1, offspring2 - SBX children (clamped to [yl,yu])
%
% Notes:
%   - Distribution index mu (a.k.a. eta_c) controls spread (larger = closer to parents).
%   - If crossover not applied, parents are returned unchanged.

    mu = 20;                                  % SBX distribution index
    offspring1 = parent1;                     % default: copy parents
    offspring2 = parent2;

    % If crossover not triggered, return clones
    if rand(1) > CrossOverRate || isempty(parent1)
        return
    end

    yl = Limits(1);                           % lower bound
    yu = Limits(2);                           % upper bound

    % Preallocate children
    child1 = zeros(size(parent1));
    child2 = zeros(size(parent2));

    % Gene-wise SBX
    for j = 1:length(parent2)
        par1 = parent1(j);
        par2 = parent2(j);

        if rand(1) <= 0.5                      % with prob 0.5 perform SBX on this gene
            if abs(par1 - par2) > 1e-6         % only if parents differ
                % Sort endpoints (y1 <= y2)
                if par2 > par1
                    y1 = par1; y2 = par2;
                else
                    y1 = par2; y2 = par1;
                end

                % Compute beta (distance to bounds), then alpha
                if (y1 - yl) > (yu - y2)
                    beta = 1 + 2*(yu - y2)/(y2 - y1);
                else
                    beta = 1 + 2*(y1 - yl)/(y2 - y1);
                end
                beta  = 1 / beta;
                alpha = 2 - beta^(mu + 1);

                % Sample spread factor betaq
                r = rand(1);
                if r <= 1/alpha
                    betaq = (r*alpha)^(1/(mu + 1));
                else
                    betaq = (1/(2 - r*alpha))^(1/(mu + 1));
                end

                % Create symmetric children around mid-point
                child1(j) = 0.5*((y1 + y2) - betaq*(y2 - y1));
                child2(j) = 0.5*((y1 + y2) + betaq*(y2 - y1));
            else
                % Parents equal: children equal to parents
                y1 = par1; y2 = par2;
                child1(j) = 0.5*((y1 + y2) - (y2 - y1));
                child2(j) = 0.5*((y1 + y2) + (y2 - y1));
            end

            % Clamp to bounds
            if child1(j) < yl, child1(j) = yl; end
            if child2(j) < yl, child2(j) = yl; end
            if child1(j) > yu, child1(j) = yu; end
            if child2(j) > yu, child2(j) = yu; end
        else
            % No per-gene crossover: copy parents on this locus
            child1(j) = par1;
            child2(j) = par2;
        end

        % Numerical safety (should not trigger under normal conditions)
        if ~isreal(child1(j)) || ~isreal(child2(j))
            % Fall back to uniform within bounds for any invalid gene
            if ~isreal(child1(j)), child1(j) = rand(1)*(yu - yl) + yl; end
            if ~isreal(child2(j)), child2(j) = rand(1)*(yu - yl) + yl; end
        end
    end

    offspring1 = child1;                       % assign results
    offspring2 = child2;
end
