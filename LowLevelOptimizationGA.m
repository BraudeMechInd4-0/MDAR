
function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationGA(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1)
%LOWLEVELOPTIMIZATIONGA Solve transfer using a simple GA; return best metrics.
% Inputs:
%   rc, vc           - current chaser state (position/velocity)
%   T                - elapsed mission time
%   RelativeEpoch    - target's relative epoch offset [s]
%   model            - target propagation model handle/data
%   rt1, vt1         - target state at reference epoch
%   T1               - planned wait before maneuver [s]
%   mc, Ac, J, Cd, Re, mu1 - physical/environment parameters
% Outputs:
%   DeltaV           - best total ΔV found
%   tm               - transfer time [s]
%   de, dm           - discrete environment/drag mode flags
%   nrev             - number of revolutions used

    nVar = 2;                % number of decision variables (time + discrete mode bits)
    % VarSize=[1 nVar];      % (unused in current implementation)
    VarMin = 10e-5;          % lower bound for real-valued variable (not enforced below)
    VarMax = 1;              % upper bound for real-valued variable

    % Pack parameters for model evaluation
    Params.rc = rc; Params.vc = vc; Params.T = T;
    Params.RelativeEpoch = RelativeEpoch;
    Params.rt1 = rt1; Params.vt1 = vt1; Params.T1 = T1;
    Params.Ac = Ac; Params.mc = mc; Params.J = J;
    Params.mu1 = mu1; Params.Cd = Cd; Params.Re = Re;

    % Maximum number of generations
    MaxIt = 20;

    % Population size (odd value)
    PopSize = (4 + round(3*log(nVar))) * 2 + 1;

    % --- Initial population ---
    ParentsX = rand(PopSize,1);          % real-valued gene in [0,1]
    ParentsI = randi([0,1],PopSize,4);   % 4 binary genes (map to 1..16 settings)

    % Evaluate initial population
    [CostParent, IDX] = EvaluateModel(ParentsX, ParentsI, Params, model);

    % Track best solution so far
    [DeltaV, I] = min(CostParent);
    Xbest   = ParentsX(I);
    IDXBest = IDX(I);

    % --- GA loop ---
    for gen = 1:MaxIt
        % Selection (binary tournament)
        [MTX, MTI] = SelectionLow(ParentsX, ParentsI, CostParent, PopSize);

        % Reproduction (crossover + mutation)
        [OffspringX, OffspringI] = ReproductionLow(MTX, MTI, PopSize);

        % Evaluate offspring
        [CostOffspring, IDX] = EvaluateModel(OffspringX, OffspringI, Params, model);

        % Update global best
        if min(CostOffspring) < DeltaV
            [DeltaV, I] = min(CostOffspring);
            Xbest   = OffspringX(I);
            IDXBest = IDX(I);
        end

        % Survivor selection (elitist replacement of worst with best)
        [ParentsX, ParentsI, CostParent] = updatePop(ParentsX, ParentsI, CostParent, ...
                                                     OffspringX, OffspringI, CostOffspring, PopSize);
    end

    % Discrete tables (drag mode, environment, revolutions)
    DM   = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
    DE   = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
    Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];

    % Map best normalized time to seconds, and choose discrete flags
    tm   = Xbest * 3 * 60 * 60;         % scale to seconds
    de   = DE(IDXBest);                 % environment flag
    dm   = DM(IDXBest);                 % drag mode flag
    nrev = Nrev(IDXBest);               % number of revolutions
end

% --- Cost evaluator: returns ΔV and chosen setting index ---
function [Y_pop,IDX] = EvaluateModel(X, Idx, Params, model)
    T  = Params.T;             RelativeEpoch = Params.RelativeEpoch;
    t1 = Params.T1;            rt1 = Params.rt1; vt1 = Params.vt1;
    rc = Params.rc;            vc = Params.vc;  mu1 = Params.mu1;
    Cd = Params.Cd;            Ac = Params.Ac;  mc = Params.mc;
    Re = Params.Re;            J  = Params.J;

    % Orbital periods (used to set integrator section size)
    [at,~,~,~,~,~] = kep_elements(rt1, vt1, mu1);
    Torbt = sqrt(at^3/mu1) * 2*pi;
    [ac,~,~,~,~,~] = kep_elements(rc, vc, mu1);
    Torbc = sqrt(ac^3/mu1) * 2*pi;

    % Integrator options for Lambert solve
    options.AbsTol = 1e-7; options.RelTol = 1e-5; options.N = 16;
    delta = 8; options.Sec = min(Torbt, Torbc) / delta;

    % Outputs
    Y_pop = nan(size(X,1),1);          % total ΔV per candidate
    IDX   = nan(size(X,1),1);          % chosen discrete setting index

    TM   = X(:,1);                     % normalized transfer time
    Sett = sum(Idx .* [1 2 4 8], 2) + 1;   % 4-bit to index (1..16)

    parfor N = 1:size(X,1)
        tm = TM(N) * 3 * 60 * 60;      % scale to seconds
        I  = Sett(N);                  % table index

        % Infeasibility: exceeds 12-month budget after waits & elapsed time
        if RelativeEpoch > 2*60*60*24*30*12 - tm - t1 - T
            Y_pop(N,1) = 1000;         % heavy penalty
        else
            % Propagate target to arrival time
            xt  = propCheb(RelativeEpoch + 1e-10 + t1 + T + tm, model);
            rt2 = xt(1:3); vt2 = xt(4:6);

            % Tables for drag/environment/revs
            DM   = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
            DE   = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
            Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];

            % Lambert transfer with J6+drag dynamics
            [vc1m, vc2m, ~] = prtlambertbf(@(t,x) orbit_eq_J6_drag(t,x,mu1,Cd,Ac,mc,Re,J), ...
                                           @odeMPCI, rc, rt2, vc, DM(I), DE(I), Nrev(I), ...
                                           tm, mu1, 1e-5, options);

            if isreal(vc1m) && isreal(vc2m)
                deltav1   = norm(vc1m - vc);   % depart ΔV
                deltav2   = norm(vc2m - vt2);  % arrive ΔV
                deltavtot = deltav1 + deltav2;
            else
                deltavtot = 1000;              % penalty if solver fails
            end

            Y_pop(N,1) = deltavtot;           % candidate cost
            IDX(N,1)   = I;                   % chosen setting index
        end
    end
end

% --- Binary tournament selection ---
function [MTX,MTI] = SelectionLow(ParentsX,ParentsI,CostParent,PopSize)
    MTX = nan(PopSize,1);
    MTI = nan(PopSize,4);
    for i = 1:PopSize
        S = randi(PopSize, 2, 1);             % pick two competitors
        % NOTE: compares their costs; winner has lower cost
        if CostParent(S(1)) < CostParent(S(2))
            MTX(i,:) = ParentsX(S(1),:);
            MTI(i,:) = ParentsI(S(1),:);
        else
            MTX(i,:) = ParentsX(S(2),:);
            MTI(i,:) = ParentsI(S(2),:);
        end
    end
end

% --- Crossover (real blend + one-point on 4-bit) and mutation ---
function [OffspringX,OffspringI] = ReproductionLow(MTX,MTI,PopSize)
    Offspring1X = MTX;                         % start from selected parents
    Offspring1I = MTI;

    %% Crossover
    pc = 0.8;                                  % crossover probability
    for i = 1:2:PopSize-1                       % pairwise (i,i+1)
        if rand() < pc
            alpha = rand();                    % blend factor
            % Real-valued blend crossover
            Offspring1X(i,:)   = MTX(i,:)*alpha + MTX(i+1,:)*(1-alpha);
            Offspring1X(i+1,:) = MTX(i+1,:)*alpha + MTX(i,:)*(1-alpha);
            % One-point crossover on 4 binary genes
            point = randi(length(Offspring1I(i,:)) - 1);
            Offspring1I(i,:)   = [MTI(i,1:point),   MTI(i+1,point+1:end)];
            Offspring1I(i+1,:) = [MTI(i+1,1:point), MTI(i,point+1:end)];
        end
    end

    % Mutation
    pm = 0.1;                                  % mutation probability
    OffspringX = Offspring1X;
    OffspringI = Offspring1I;
    for i = 1:PopSize
        if rand < pm
            OffspringX(i,:) = rand(size(Offspring1X(i,:)));     % re-sample real gene
            point = randi(length(OffspringI(i,:)));             % flip one bit
            OffspringI(i,point) = 1 - OffspringI(i,point);
        end
    end
end

% --- Elitist survivor selection: keep best PopSize ---
function [NextParentsX,NextParentsI,NextParentCost] = updatePop(ParentsX,ParentsI,CostParent,OffspringX,OffspringI,CostOffspring,PopSize)
    PopAllX = [ParentsX; OffspringX];          % merge parents + offspring
    PopAllI = [ParentsI; OffspringI];
    CostAll = [CostParent; CostOffspring];
    [~, idx] = sort(CostAll);                  % sort ascending by cost
    NextParentsX    = PopAllX(idx(1:PopSize),:);
    NextParentsI    = PopAllI(idx(1:PopSize),:);
    NextParentCost  = CostAll(idx(1:PopSize));
end
