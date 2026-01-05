
function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationRandom(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1)
%LOWLEVELOPTIMIZATIONRANDOM Sample random candidates; pick the best ΔV.
% Inputs:
%   rc, vc           - current chaser state (position/velocity)
%   T                - elapsed mission time
%   RelativeEpoch    - target's relative epoch offset [s]
%   model            - target propagation model
%   rt1, vt1         - target state at reference epoch
%   T1               - planned wait before maneuver [s]
%   mc, Ac, J, Cd, Re, mu1 - physical/environment parameters
% Outputs:
%   DeltaV           - best total ΔV found among random samples
%   tm               - corresponding transfer time [s]
%   de, dm           - discrete environment/drag mode flags
%   nrev             - number of revolutions used

    nVar  = 2;                   % number of decision variables (time + 4-bit mode)
    VarMin = 1e-5;               % lower bound for normalized time (unused in sampling)
    VarMax = 1;                  % upper bound for normalized time

    % Pack parameters for evaluator
    Params.rc = rc; Params.vc = vc; Params.T = T;
    Params.RelativeEpoch = RelativeEpoch;
    Params.rt1 = rt1; Params.vt1 = vt1; Params.T1 = T1;
    Params.Ac = Ac; Params.mc = mc; Params.J = J;
    Params.mu1 = mu1; Params.Cd = Cd; Params.Re = Re;

    MaxIt  = 20;                 % number of random batches
    lambda = (4 + round(3*log(nVar))) * 2;   % batch size per iteration

    % Total number of random solutions to evaluate
    Nsolutions = MaxIt * (lambda + 1) + 1;

    % --- Random population: times in [0,1], 4 bits for discrete setting ---
    X   = rand(Nsolutions, 1);               % normalized transfer times
    Idx = randi([0,1], Nsolutions, 4);       % 4 binary genes (map to 1..16)

    % Evaluate all candidates
    [Cost, IDX] = EvaluateModel(X, Idx, Params, model);

    % Pick best by minimum cost (ΔV)
    [DeltaV, I] = min(Cost);

    % Discrete tables (drag mode, environment, revolutions)
    DM   = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
    DE   = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
    Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];

    % Map best normalized time to seconds, and choose discrete flags
    tm   = X(I) * 3 * 60 * 60;              % scale (note: 3h per unit)
    de   = DE(IDX(I));
    dm   = DM(IDX(I));
    nrev = Nrev(IDX(I));
end

% --- Cost evaluator: returns ΔV and chosen setting index ---
function [Y_pop,IDX] = EvaluateModel(X,Idx,Params,model)

    % Unpack parameters
    T  = Params.T;             RelativeEpoch = Params.RelativeEpoch;
    t1 = Params.T1;            rt1 = Params.rt1; vt1 = Params.vt1;
    rc = Params.rc;            vc = Params.vc;  mu1 = Params.mu1;
    Cd = Params.Cd;            Ac = Params.Ac;  mc = Params.mc;
    Re = Params.Re;            J  = Params.J;

    % Orbital periods (for integrator section size)
    [at,~,~,~,~,~] = kep_elements(rt1, vt1, mu1);
    Torbt = sqrt(at^3/mu1) * 2*pi;
    [ac,~,~,~,~,~] = kep_elements(rc, vc, mu1);
    Torbc = sqrt(ac^3/mu1) * 2*pi;

    % Integrator options for Lambert solve
    options.AbsTol = 1e-7; options.RelTol = 1e-5; options.N = 16;
    delta = 8; options.Sec = min(Torbt, Torbc) / delta;

    % Outputs
    Y_pop = nan(size(X,1),1);              % total ΔV per candidate
    IDX   = nan(size(X,1),1);              % chosen discrete setting index

    TM   = X(:,1);                         % normalized transfer time
    Sett = sum(Idx .* [1 2 4 8], 2) + 1;   % 4-bit to index (1..16)

    % Evaluate each candidate (parallel-friendly)
    parfor N = 1:size(X,1)
        tm = TM(N) * 3 * 60 * 60;          % scale to seconds (note: 3h per unit)
        I  = Sett(N);                       % table index

        % Infeasibility: exceeds 12-month budget after waits & elapsed time
        if RelativeEpoch > 2*60*60*24*30*12 - tm - t1 - T
            Y_pop(N,1) = 1000;             % heavy penalty
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
