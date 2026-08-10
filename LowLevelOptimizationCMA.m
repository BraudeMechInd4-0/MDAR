
function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationCMA(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1)
%LOWLEVELOPTIMIZATIONCMA Solve transfer using CMA-ES and return best metrics.
% Inputs:
%   rc, vc           - current chaser state (pos/vel)
%   T                - elapsed mission time
%   RelativeEpoch    - target's relative epoch offset [s]
%   model            - target propagation model
%   rt1, vt1         - target state at reference epoch
%   T1               - planned wait before maneuver [s]
%   mc, Ac, J, Cd, Re, mu1 - physical/environment parameters
% Outputs:
%   DeltaV           - total ΔV of the selected solution
%   tm               - transfer time [s]
%   de, dm           - discrete environment/drag mode flags
%   nrev             - number of revolutions used

    nVar   = 2;                 % number of decision variables (tm, setting)
    VarSize = [1 nVar];         % decision vector shape
    VarMin  = 0;                % lower bound (normalized)
    VarMax  = 1;                % upper bound (normalized)

    % Pack parameters for evaluation
    Params.rc = rc; Params.vc = vc; Params.T = T;
    Params.RelativeEpoch = RelativeEpoch;
    Params.rt1 = rt1; Params.vt1 = vt1; Params.T1 = T1;
    Params.Ac = Ac; Params.mc = mc; Params.J = J;
    Params.mu1 = mu1; Params.Cd = Cd; Params.Re = Re;

    %% CMA-ES settings
    MaxIt  = 10;                                                % max iterations
    lambda = 2*(4 + round(3*log(nVar))) * 2;                    % offspring size
    mu     = round(lambda/2);                                    % parent count
    w      = log(mu + 0.5) - log(1:mu); w = w/sum(w);           % parent weights
    mu_eff = 1/sum(w.^2);                                       % effective μ

    % Step-size control
    sigma0 = 0.3*(VarMax - VarMin);
    cs = (mu_eff + 2)/(nVar + mu_eff + 5);
    ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar + 1)) - 1, 0);
    ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));

    % Covariance update parameters
    cc  = (4 + mu_eff/nVar)/(4 + nVar + 2*mu_eff/nVar);
    c1  = 2/((nVar + 1.3)^2 + mu_eff);
    alpha_mu = 2;
    cmu = min(1 - c1, alpha_mu*(mu_eff - 2 + 1/mu_eff)/((nVar + 2)^2 + alpha_mu*mu_eff/2));
    hth = (1.4 + 2/(nVar + 1)) * ENN;

    % Initialize paths
    ps    = cell(MaxIt,1); pc = cell(MaxIt,1);
    C     = cell(MaxIt,1); sigma = cell(MaxIt,1);
    ps{1} = zeros(VarSize); pc{1} = zeros(VarSize);
    C{1}  = eye(nVar);     sigma{1} = sigma0;

    % Individual template
    empty_individual.Position = [];
    empty_individual.Step     = [];
    empty_individual.Cost     = [];

    % Initialize mean (M) and best
    M = repmat(empty_individual, MaxIt, 1);
    M(1).Position = unifrnd(VarMin, VarMax, VarSize);           % random start in [0,1]
    M(1).Step     = zeros(VarSize);
    [M(1).Cost, M(1).IDX] = EvaluateModel(M(1).Position, Params, model);
    BestSol = M(1);                                             % best-so-far
    BestCost = zeros(MaxIt,1);                                  % log best cost

    %% CMA-ES main loop
    for g = 1:MaxIt
        % --- Sample offspring around current mean ---
        pop = repmat(empty_individual, lambda, 1);
        X   = nan(lambda, nVar);                                % store positions for batch eval
        for i = 1:lambda
            pop(i).Step     = mvnrnd(zeros(VarSize), C{g});     % sample step ~ N(0, C)
            pop(i).Position = M(g).Position + sigma{g} * pop(i).Step; % new candidate
            pop(i).Position = max([pop(i).Position; 1e-10*ones(size(pop(i).Position))], [], 1); % clamp low
            pop(i).Position = min([pop(i).Position; ones(size(pop(i).Position))], [], 1);       % clamp high
            X(i,:) = pop(i).Position;
        end

        % --- Evaluate all offspring at once ---
        [CurrentCosts, CurrentIDX] = EvaluateModel(X, Params, model);
        for i = 1:lambda
            pop(i).Cost = CurrentCosts(i);
            pop(i).IDX  = CurrentIDX(i);
            if pop(i).Cost < BestSol.Cost                        % update global best
                BestSol = pop(i);
            end
        end

        % --- Sort by cost (ascending) and save best cost ---
        Costs = [pop.Cost];
        [~, SortOrder] = sort(Costs);
        pop      = pop(SortOrder);
        BestCost(g) = BestSol.Cost;

        % --- Stop at last iteration ---
        if g == MaxIt
            break
        end

        % --- Update mean step and position ---
        M(g+1).Step = 0;
        for j = 1:mu
            M(g+1).Step = M(g+1).Step + w(j) * pop(j).Step;
        end
        M(g+1).Position = M(g).Position + sigma{g} * M(g+1).Step;
        M(g+1).Position = max([M(g+1).Position; 1e-20*ones(size(M(g+1).Position))], [], 1);
        M(g+1).Position = min([M(g+1).Position; ones(size(M(g+1).Position))], [], 1);
        [M(g+1).Cost, M(g+1).IDX] = EvaluateModel(M(g+1).Position, Params, model);
        if M(g+1).Cost < BestSol.Cost
            BestSol = M(g+1);
        end

        % --- Update step-size (σ) ---
        ps{g+1}    = (1 - cs)*ps{g} + sqrt(cs*(2 - cs)*mu_eff)*M(g+1).Step/chol(C{g})';
        sigma{g+1} = sigma{g} * exp(cs/ds*(norm(ps{g+1})/ENN - 1))^0.3;

        % --- Update covariance (C) ---
        if norm(ps{g+1})/sqrt(1 - (1 - cs)^(2*(g+1))) < hth
            hs = 1;
        else
            hs = 0;
        end
        delta = (1 - hs)*cc*(2 - cc);
        pc{g+1} = (1 - cc)*pc{g} + hs*sqrt(cc*(2 - cc)*mu_eff)*M(g+1).Step;
        C{g+1}  = (1 - c1 - cmu)*C{g} + c1*(pc{g+1}'*pc{g+1} + delta*C{g});
        for j = 1:mu
            C{g+1} = C{g+1} + cmu*w(j)*pop(j).Step'*pop(j).Step;
        end

        % --- Ensure C is positive semidefinite ---
        [V, E] = eig(C{g+1});
        if any(diag(E) < 0)
            E = max(E, 0);
            C{g+1} = V*E/V;
        end
    end

    % --- Map best solution to outputs ---
    DeltaV = BestSol.Cost;                                      % total ΔV
    DM  = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S']; % drag mode table
    DE  = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L']; % environment table
    Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];                                % revolutions table

    tm  = BestSol.Position(1) * 3 * 60 * 60;                   % scale normalized time -> seconds
    de  = DE(BestSol.IDX);                                     % selected environment flag
    dm  = DM(BestSol.IDX);                                     % selected drag mode
    nrev = Nrev(BestSol.IDX);                                  % selected revolutions
end

% --- Batch evaluator used by CMA-ES ---
function [Y_pop,IDX] = EvaluateModel(X, Params, model)
%EVALUATEMODEL Evaluate ΔV for candidates X and return cost + chosen setting.
% Inputs:
%   X        - [N x 2] (tm_norm, setting_norm)
%   Params   - packed parameters (states, env, timing)
%   model    - target propagation model
% Outputs:
%   Y_pop    - ΔV per candidate (penalized if infeasible)
%   IDX      - chosen discrete setting index per candidate (1..16)

    % Unpack parameters
    T  = Params.T;            RelativeEpoch = Params.RelativeEpoch;
    t1 = Params.T1;           rt1 = Params.rt1; vt1 = Params.vt1;
    rc = Params.rc;           vc = Params.vc;  mu1 = Params.mu1;
    Cd = Params.Cd;           Ac = Params.Ac;  mc = Params.mc;
    Re = Params.Re;           J  = Params.J;

    % Orbital periods for step sizing
    [at,~,~,~,~,~] = kep_elements(rt1, vt1, mu1);
    Torbt = sqrt(at^3/mu1)*2*pi;
    [ac,~,~,~,~,~] = kep_elements(rc, vc, mu1);
    Torbc = sqrt(ac^3/mu1)*2*pi;

    % Integrator/options for Lambert + propagation
    options.AbsTol = 1e-7; options.RelTol = 1e-5; options.N = 16;
    delta = 8; options.Sec = min(Torbt, Torbc)/delta;

    % Prepare outputs
    Y_pop = nan(size(X,1),1);
    IDX   = nan(size(X,1),1);

    TM   = X(:,1);                  % normalized transfer time
    Sett = round(X(:,2)*15) + 1;    % discrete setting index 1..16

    % Evaluate each candidate (parallel friendly)
    parfor N = 1:size(X,1)
        tm = TM(N) * 3 * 60 * 60;   % scale to seconds
        I  = Sett(N);               % table index

        % Infeasibility: excessive time vs. 12 months budget
        if RelativeEpoch > 2*60*60*24*30*12 - tm - t1 - T
            Y_pop(N,1) = 1000;      % heavy penalty
        else
            % Target propagation to arrival time
            xt  = propCheb(RelativeEpoch + 1e-10 + t1 + T + tm, model);
            rt2 = xt(1:3); vt2 = xt(4:6);

            % Predefined tables (drag/environment/revs)
            DM  = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
            DE  = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
            Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];

            % Lambert solve with J6+drag dynamics
            [vc1m, vc2m] = prtlambertT(@(t,x) orbit_eq_J6_drag(t,x,mu1,Cd,Ac,mc,Re,J), ...
                                       @odeMPCI, rc, rt2, vc, DM(I), DE(I), Nrev(I), ...
                                       tm, 1, options);

            if isreal(vc1m) && isreal(vc2m)
                deltav1   = norm(vc1m - vc);  % depart ΔV
                deltav2   = norm(vc2m - vt2); % arrive ΔV
                deltavtot = deltav1 + deltav2;
            else
                deltavtot = 1000;             % penalty if solver fails
            end

            Y_pop(N,1) = deltavtot;           % cost for candidate N
            IDX(N,1)   = I;                   % chosen setting index
        end
    end
end
