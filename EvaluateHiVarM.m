function PopulationAfter = EvaluateHiVarM(Population,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLeverType,Solver)
%EVALUATEHIVARM  High-level evaluation with a mass-varying chaser (multi-objective).
%   PopulationAfter = EVALUATEHIVARM(Population, model, m, RelativeEpoch, mc, Ac, J,
%   Cd, Re, mu1, wM, wDV, LowLeverType, Solver) evaluates each removal plan by running
%   the lower-level solver on every leg while decreasing the chaser mass as propellant
%   is expended, recording removed mass, per-leg DeltaV, and timing. This is the
%   evaluator used by the NSGA-II main and the weighted-sum baseline.
%
%   Inputs:
%     Population    - struct array with fields .Order, .WaitUntilManuver
%     model         - per-debris Chebyshev propagation models
%     m             - debris mass vector (indexed by .Order)
%     RelativeEpoch - per-debris relative epoch offsets [s]
%     mc, Ac        - initial chaser mass [kg] and cross-sectional area [km^2]
%     J             - zonal harmonic coefficients [J1..J6]
%     Cd, Re, mu1   - drag coefficient, Earth radius [km], mu [km^3/s^2]
%     wM, wDV       - weights for the scalarized .fit field (unused by NSGA-II)
%     LowLeverType  - lower-level optimizer: "Rand" | "GA" (any other value uses
%                     the one-variable optimizer)
%     Solver        - Lambert-solver select passed to LowLevelOptimizationRandom
%   Outputs:
%     PopulationAfter - input augmented with .DV, .DE, .TM, .DM, .NREV,
%                       .Nremoved, .Mremoved, .WaitUntilManuver, .fit
mc0 = mc;
DeltaVmax = 100;
PopulationAfter = Population;
for i = 1 : length(Population)

    Order = Population(i).Order;

    T1 = Population(i).WaitUntilManuver;
    if length(Order) == length(unique(Order))
        xc =  propCheb(RelativeEpoch(Order(1))+1e-10+T1(1),model(Order(1)));
        rc = xc(1:3);
        vc = xc(4:6);
        Mremoved = m(Order(1));
        Nremoved = 1;
        T = 0;
        DV = nan(length(Order) - 1,1);
        DE = nan(length(Order) - 1,1);
        TM = zeros(length(Order) - 1,1);
        DM = nan(length(Order) - 1,1);
        NREV = nan(length(Order) - 1,1);
        for ManIdx = 1:length(Order) - 1
            idxT = Order(ManIdx + 1);
            T = T + T1(ManIdx);
            xt = propCheb(RelativeEpoch(idxT)+1e-10+T,model(idxT));
            rt1 = xt(1:3);
            vt1 = xt(4:6);
            if strcmp(LowLeverType,"Rand")
                [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationRandom(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,0,mc,Ac,J,Cd,Re,mu1,Solver);
            elseif strcmp(LowLeverType,"GA")
                [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationGA(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,0,mc,Ac,J,Cd,Re,mu1);
            else
                [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationOneVar(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,0,mc,Ac,J,Cd,Re,mu1);
            end

            % Varing M
            mc = min(mc0/2,mc - DeltaV/DeltaVmax*mc0/2);
            T = T + tm;
            if ManIdx<length(Order) - 1
            xc = propCheb(RelativeEpoch(idxT)+1e-10+T+T1(ManIdx+1),model(idxT));
            rc = xc(1:3);
            vc = xc(4:6);
            end
            DV(ManIdx,1) = DeltaV;
            Nremoved = Nremoved + 1;
            Mremoved = Mremoved + m(idxT);
            if DeltaV >= 1000 || isnan(DeltaV)

                DV(ManIdx:end) = 1000;
                break
            end

            
            DE(ManIdx,1) = de;
            TM(ManIdx,1) = tm;
            DM(ManIdx,1) = dm;
            NREV(ManIdx,1) = nrev;
        end

        if  Nremoved < length(Order)
            Penalty = 1000*(length(Order)-Nremoved);
        else
            Penalty = 0;
        end
        PopulationAfter(i).DV = DV + Penalty ;
        PopulationAfter(i).Nremoved = Nremoved;
        PopulationAfter(i).Mremoved = Mremoved-Penalty^2;
        PopulationAfter(i).DE = DE;
        PopulationAfter(i).TM = TM + Penalty;
        if size(PopulationAfter(i).TM,2) > 1
                PopulationAfter(i).TM = [PopulationAfter(i).TM]';
        end
        PopulationAfter(i).DM = DM;
        PopulationAfter(i).NREV = NREV;
        PopulationAfter(i).Order = Order;
        PopulationAfter(i).WaitUntilManuver = Population(i).WaitUntilManuver + Penalty;
        if size(PopulationAfter(i).WaitUntilManuver,2) > 1
            PopulationAfter(i).WaitUntilManuver = [PopulationAfter(i).WaitUntilManuver]';
        end
        PopulationAfter(i).fit = wM*Mremoved + wDV*sum(DV);



    else

        Penalty = 1000*(length(Order));
        PopulationAfter(i).DV = zeros(length(Order) - 1,1) + Penalty^2;
        PopulationAfter(i).Nremoved = 0;
        PopulationAfter(i).Mremoved = 0;
        PopulationAfter(i).DE = [];
        PopulationAfter(i).TM = zeros(length(Order) - 1,1) + Penalty^2;
        if size(PopulationAfter(i).TM,2) > 1
                PopulationAfter(i).TM = [PopulationAfter(i).TM]';
        end
        PopulationAfter(i).DM = [];
        PopulationAfter(i).NREV = [];
        PopulationAfter(i).Order = Order;
        PopulationAfter(i).WaitUntilManuver = Population(i).WaitUntilManuver + Penalty^2;
        if size(PopulationAfter(i).WaitUntilManuver,2) > 1
            PopulationAfter(i).WaitUntilManuver = [PopulationAfter(i).WaitUntilManuver]';
        end
        PopulationAfter(i).fit = wM*0 + wDV*1000*(length(Order)-1);
    end

end

