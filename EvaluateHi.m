
function PopulationAfter = EvaluateHi(Population, model, m, RelativeEpoch, mc, Ac, J, Cd, Re, mu1, wM, wDV, LowLeverType)
%EVALUATEHI Evaluate high-level individual by chaining low-level optimizations.
% Inputs:
%   Population     - struct array with fields .Order, .WaitUntilManuver
%   model          - per-debris propagation/chebyshev model structs
%   m              - debris mass vector (aligned with indices in Order)
%   RelativeEpoch  - per-debris relative epoch offsets [s]
%   mc, Ac, J, Cd, Re, mu1 - physical/propagation parameters
%   wM, wDV        - fitness weights (mass removed vs. total ΔV)
%   LowLeverType   - 'Rand' | 'GA' | 'CMA' (choose low-level optimizer)
% Output:
%   PopulationAfter - same struct array with added metrics:
%                     .DV, .DE, .TM, .DM, .NREV, .Nremoved, .Mremoved, .fit

    PopulationAfter = Population;                             % start from input

    for i = 1:length(Population)
        Order = Population(i).Order;                          % sequence of debris indices
        T1    = Population(i).WaitUntilManuver;               % waits before each maneuver

        %--- Initialize current chaser state at first target epoch ---
        xc = propCheb(RelativeEpoch(Order(1)) + 1e-10, model(Order(1)));
        rc = xc(1:3);                                         % position at first target
        vc = xc(4:6);                                         % velocity at first target

        %--- Initialize accumulation metrics ---
        Mremoved = m(Order(1));                               % mass removed (start with first)
        Nremoved = 1;                                         % number removed so far
        T        = 0;                                         % elapsed mission time

        % Preallocate per-maneuver logs (length = number of transfers)
        DV   = nan(length(Order) - 1, 1);                     % ΔV per transfer
        DE   = nan(length(Order) - 1, 1);                     % energy change (from low-level)
        TM   = nan(length(Order) - 1, 1);                     % transfer time
        DM   = nan(length(Order) - 1, 1);                     % mass change (if applicable)
        NREV = nan(length(Order) - 1, 1);                     % revs used (if applicable)

        %--- Iterate over remaining targets in the plan ---
        for ManIdx = 1:length(Order) - 1
            idxT = Order(ManIdx + 1);                         % next target index

            % Propagate target to current mission time (T) and its relative epoch
            xt  = propCheb(RelativeEpoch(idxT) + 1e-10 + T, model(idxT));
            rt1 = xt(1:3);
            vt1 = xt(4:6);

            % Select low-level optimizer based on requested type
            if strcmp(LowLeverType, "Rand")
                [DeltaV, tm, de, dm, nrev] = LowLevelOptimizationRandom( ...
                    rc, vc, T, RelativeEpoch(idxT), model(idxT), rt1, vt1, T1(ManIdx), mc, Ac, J, Cd, Re, mu1);
            elseif strcmp(LowLeverType, "GA")
                [DeltaV, tm, de, dm, nrev] = LowLevelOptimizationGA( ...
                    rc, vc, T, RelativeEpoch(idxT), model(idxT), rt1, vt1, T1(ManIdx), mc, Ac, J, Cd, Re, mu1);
            else
                [DeltaV, tm, de, dm, nrev] = LowLevelOptimizationCMA( ...
                    rc, vc, T, RelativeEpoch(idxT), model(idxT), rt1, vt1, T1(ManIdx), mc, Ac, J, Cd, Re, mu1);
            end

            % Update chaser state at arrival (wait + transfer time)
            xc = propCheb(RelativeEpoch(idxT) + 1e-10 + T + tm, model(idxT));
            rc = xc(1:3);
            vc = xc(4:6);

            % Log metrics and update counters
            DV(ManIdx,1) = DeltaV;
            Nremoved      = Nremoved + 1;
            Mremoved      = Mremoved + m(idxT);

            %--- Failure/penalty guard: abort if ΔV invalid/too large ---
            if DeltaV >= 1000 || isnan(DeltaV)
                Nremoved = 0;
                Mremoved = 0;
                DV       = ones(size(DV)) * 1000;            % heavy penalty on fitness
                break
            end

            % Advance mission time: scheduled wait + transfer time
            T = T + T1(ManIdx) + tm;

            % Store additional low-level outputs
            DE(ManIdx,1)  = de;
            TM(ManIdx,1)  = tm;
            DM(ManIdx,1)  = dm;
            NREV(ManIdx,1)= nrev;
        end

        %--- Attach results back to individual ---
        PopulationAfter(i).DV               = DV;
        PopulationAfter(i).Nremoved         = Nremoved;
        PopulationAfter(i).Mremoved         = Mremoved;
        PopulationAfter(i).DE               = DE;
        PopulationAfter(i).TM               = TM;
        PopulationAfter(i).DM               = DM;
        PopulationAfter(i).NREV             = NREV;
        PopulationAfter(i).Order            = Order;
        PopulationAfter(i).WaitUntilManuver = T1;

        % Fitness: reward removed mass, penalize total ΔV
        PopulationAfter(i).fit = wM * Mremoved - wDV * sum(DV);
    end
end
