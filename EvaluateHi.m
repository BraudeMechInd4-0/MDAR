function PopulationAfter = EvaluateHi(Population,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLeverType)
%EVALUATEHI  Single-objective high-level evaluation via the lower-level solver.
%   PopulationAfter = EVALUATEHI(Population, model, m, RelativeEpoch, mc, Ac, J, Cd,
%   Re, mu1, wM, wDV, LowLeverType) evaluates each removal plan by running the
%   lower-level trajectory optimizer on every leg, accumulating removed mass and
%   total DeltaV, and setting the weighted-sum fitness .fit = wM*Mremoved - wDV*sum(DV).
%
%   Inputs:
%     Population    - struct array with fields .Order, .WaitUntilManuver
%     model         - per-debris Chebyshev propagation models
%     m             - debris mass vector (indexed by .Order)
%     RelativeEpoch - per-debris relative epoch offsets [s]
%     mc, Ac        - chaser mass [kg] and cross-sectional area [km^2]
%     J             - zonal harmonic coefficients [J1..J6]
%     Cd, Re, mu1   - drag coefficient, Earth radius [km], mu [km^3/s^2]
%     wM, wDV       - fitness weights (removed mass vs. total DeltaV)
%     LowLeverType  - lower-level optimizer: "Rand" | "GA" | "CMA"
%                     (any other value uses LowLevelOptimizationOneVar)
%   Outputs:
%     PopulationAfter - input augmented with per-plan metrics
%                       .DV, .DE, .TM, .DM, .NREV, .Nremoved, .Mremoved, .fit

PopulationAfter = Population;
for i = 1 : length(Population)
    Order = Population(i).Order;
    T1 = Population(i).WaitUntilManuver;
    xc =  propCheb(RelativeEpoch(Order(1))+1e-10,model(Order(1)));
    rc = xc(1:3);
    vc = xc(4:6);
    Mremoved = m(Order(1));
    Nremoved = 1;
    T = 0;
    DV = nan(length(Order) - 1,1);
    DE = nan(length(Order) - 1,1);
    TM = nan(length(Order) - 1,1);
    DM = nan(length(Order) - 1,1);
    NREV = nan(length(Order) - 1,1);
    for ManIdx = 1:length(Order) - 1
        idxT = Order(ManIdx + 1);
        xt = propCheb(RelativeEpoch(idxT)+1e-10+T,model(idxT));
        rt1 = xt(1:3);
        vt1 = xt(4:6);
        if strcmp(LowLeverType,"Rand")
            [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationRandom(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,T1(ManIdx),mc,Ac,J,Cd,Re,mu1,1);
        elseif strcmp(LowLeverType,"GA")
            [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationGA(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,T1(ManIdx),mc,Ac,J,Cd,Re,mu1);
        elseif strcmp(LowLeverType,"CMA")
            [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationCMA(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,T1(ManIdx),mc,Ac,J,Cd,Re,mu1);
        else
            [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationOneVar(rc,vc,T,RelativeEpoch(idxT),model(idxT),rt1,vt1,T1(ManIdx),mc,Ac,J,Cd,Re,mu1);
        end

        xc = propCheb(RelativeEpoch(idxT)+1e-10+T+tm,model(idxT));
        rc = xc(1:3);
        vc = xc(4:6);
        DV(ManIdx,1) = DeltaV;
        Nremoved = Nremoved + 1;
        Mremoved = Mremoved + m(idxT);
        if DeltaV >= 1000 || isnan(DeltaV) 
            Nremoved = 0;
            Mremoved = 0;
            DV = ones(size(DV))*1000;
            break
        end
        
        T = T + T1(ManIdx) + tm;
        DE(ManIdx,1) = de;
        TM(ManIdx,1) = tm;
        DM(ManIdx,1) = dm;
        NREV(ManIdx,1) = nrev;
    end
    PopulationAfter(i).DV = DV;
    PopulationAfter(i).Nremoved = Nremoved;
    PopulationAfter(i).Mremoved = Mremoved;
    PopulationAfter(i).DE = DE;
    PopulationAfter(i).TM = TM;
    PopulationAfter(i).DM = DM;
    PopulationAfter(i).NREV = NREV;
    PopulationAfter(i).Order = Order;
    PopulationAfter(i).WaitUntilManuver = T1;
    PopulationAfter(i).fit = wM*Mremoved - wDV*sum(DV);

end

