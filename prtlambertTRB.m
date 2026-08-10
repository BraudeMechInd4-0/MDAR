function [v1t, v2t] = prtlambertTRB(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
%PRTLAMBERTTB  Perturbed Lambert solver — Thompson (select fwd/bwd) then Broyden.
%
% Preamble: propagate forward from r1 and backward from r2; choose the
%           direction with the smaller initial residual.
% Phase 1:  Thompson shooting in the chosen direction for N_switch iterations
%           or until normdr < T_switch.
% Phase 2:  Broyden's quasi-Newton method (forward), initialized with a
%           secant from the lambertb seed to the Thompson exit state.
%
% History always tracks the FORWARD error (normdr2), even during backward phase.
%
% Usage:
%   [v1t,v2t,hitearthvar,hitrad,integration_number,history] = prtlambertTRB( ...
%       forcemodel, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)

mu = 398600.5;
r_GEO = 42164;
r_max_allowed = 1.5 * r_GEO;
flag = 0;
if nargin < 9
    Learning_rate = 1;
end
if nargin < 10,  options = struct();  end
if ~isfield(options, 'Tolr'),    options.Tolr    = 1e-5;  end
if ~isfield(options, 'maxIter'), options.maxIter = 100;   end
if ~isfield(options, 'delta'),   options.delta   = 16;    end

Tolr    = options.Tolr;
maxIter = options.maxIter;
debug   = isfield(options, 'debug') && options.debug;

T_switch = Tolr * 1000;
N_switch = floor(maxIter * 0.05);
if isfield(options, 'T_switch'),  T_switch = options.T_switch;  end
if isfield(options, 'N_switch'),  N_switch = options.N_switch;  end

if debug
    fprintf('starting TRB\n');
    history = struct('iterations', [], 'normdr2', [], 'phase', {{}}, ...
        'direction_chosen', '', 'converged', false, 'final_error', NaN);
else
    history = struct();
end

%% ---- Preamble: lambertb + forward + backward propagations ----
[v1tnew, v2tnew] = lambertb(r1, r2, v1, dm, de, nrev, dtsec,mu);
if ~all(isreal(v1tnew)) || any(isnan(v1tnew)) || ~all(isreal(v2tnew)) || any(isnan(v2tnew))
    v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
    hitrad = NaN;  hitearthvar = NaN;
    integration_number = 0;
    return
end
if  hitearth ( 100, r1, v1tnew, r2, v2tnew, nrev )
    v1t = v1tnew + 1i;
    v2t = v2tnew + 1i;
    errorb = 'Hit earth';
    return
end

%% Forward propagation (from r1)
oe = kep_elements(r1, v1tnew, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_fwd     = options;
options_fwd.Sec = T_orb / options.delta;
[~, xoutfwd]    = odeMPCI(forcemodel, [0 dtsec], [r1 v1tnew], options_fwd);

%% Backward propagation (from r2)
oe = kep_elements(r2, v2tnew, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_bkw     = options;
options_bkw.Sec = T_orb / options.delta;
[~, xoutbkw]    = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew], options_bkw);

integration_number = 2;
normdr1 = norm(xoutbkw(1,1:3) - r1);
normdr2 = norm(xoutfwd(end,1:3) - r2);

%% Save lambertb seed for Broyden initialization
v1_A   = v1tnew;
dr2_A  = xoutfwd(end,1:3) - r2;   % unscaled forward error from lambertb

if debug
    history.iterations(end+1) = integration_number;
    history.normdr2(end+1)    = normdr2;
    history.phase{end+1}      = 'init';
    fprintf('init (TRB): normdr1=%g, normdr2=%g\n', normdr1, normdr2);
end

%% ---- Select direction ----
use_backward = isfinite(normdr1) && (normdr1 < normdr2);
if debug
    if use_backward
        history.direction_chosen = 'backward';
        fprintf('TRB: using backward\n');
    else
        history.direction_chosen = 'forward';
        fprintf('TRB: using forward\n');
    end
end

%% ====================================================================
%% BACKWARD BRANCH
%% ====================================================================
if use_backward
    prop_r1 = r1;
    deltar1 = (xoutbkw(1,1:3) - r1) * Learning_rate;
    xoutbkw_cur = xoutbkw;
    v2tnew_cur  = v2tnew;
    J_rev      = eye(3);
    v2_t_prev  = v2tnew;
    dr1_t_prev = xoutbkw(1,1:3) - r1;

    n_thompson   = 0;
    normdr1_prev = norm(deltar1);
    while norm(deltar1) > T_switch && (n_thompson < N_switch || norm(deltar1) > normdr1_prev) && integration_number < maxIter
        prop_r1 = prop_r1 - deltar1;
        v1tmp   = xoutbkw_cur(1,4:6);
        [v1tmp, v2tnew_cur] = lambertb(prop_r1, r2, v1tmp, dm, de, nrev, dtsec,mu);
        if ~all(isreal(v1tmp)) || any(isnan(v1tmp)) || ~all(isreal(v2tnew_cur)) || any(isnan(v2tnew_cur))
            if debug,  history.converged = false;  end
             v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
            hitrad = NaN;  hitearthvar = NaN;
            return
        end

        if  hitearth ( 100, r1, v1tmp, r2, v2tnew_cur, nrev )
             v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
            errorb = 'Hit earth';
            return
        end

        oe = kep_elements(r2, v2tnew_cur, mu);
        a  = oe(1);
        if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
            T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_orb = sqrt(a^3/mu) * 2*pi;
        end
        options_bkw.Sec = T_orb / options.delta;
        [~, xoutbkw_cur] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew_cur], options_bkw);
        integration_number = integration_number + 1;
        if ~all(isfinite(xoutbkw_cur(1,:))), break, end
        normdr1_prev = norm(deltar1);
        deltar1   = (xoutbkw_cur(1,1:3) - r1) * Learning_rate;
        v2_t_cur  = v2tnew_cur;
        dr1_t_cur = xoutbkw_cur(1,1:3) - r1;
        dx_r = (v2_t_cur  - v2_t_prev)';
        df_r = (dr1_t_cur - dr1_t_prev)';
        if norm(dx_r) > 0
            J_rev = J_rev + (df_r - J_rev*dx_r) * dx_r' / (dx_r'*dx_r);
        end
        v2_t_prev  = v2_t_cur;
        dr1_t_prev = dr1_t_cur;
        n_thompson = n_thompson + 1;

        if debug
            %% Track forward error (extra integration)
            v1_check = xoutbkw_cur(1,4:6);
            oe_c = kep_elements(r1, v1_check, mu);
            a_c  = oe_c(1);
            if a_c <= 0 || ~isreal(a_c) || isnan(a_c) || a_c > 2*r_max_allowed
                T_c = 2*pi*sqrt((2*r_max_allowed)^3/mu);
            else
                T_c = sqrt(a_c^3/mu) * 2*pi;
            end
            opt_c = options;  opt_c.Sec = T_c / options.delta;
            [~, xfwd_c] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_check], opt_c);
            integration_number = integration_number + 1;
            history.iterations(end+1) = integration_number;
            history.normdr2(end+1)    = norm(xfwd_c(end,1:3) - r2);
            history.phase{end+1}      = 'Thompson_bkw';
            fprintf('integration %d (T bkw): normdr1=%g, normdr2=%g\n', ...
                integration_number, norm(deltar1)/Learning_rate, norm(xfwd_c(end,1:3)-r2));
        end
    end

    %% Get v1_B and dr2_B after backward Thompson phase
    v1_B = xoutbkw_cur(1,4:6);
    oe = kep_elements(r1, v1_B, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_fwd.Sec = T_orb / options.delta;
    [~, xout_B]    = odeMPCI(forcemodel, [0 dtsec], [r1 v1_B], options_fwd);
    integration_number = integration_number + 1;
    dr2_B   = xout_B(end,1:3) - r2;
    xout    = xout_B;
    v1_cur  = v1_B;
    dr2_cur = dr2_B;

    J = -J_rev';   % ∂r2/∂v1 ≈ -(∂r1/∂v2)^T (symplectic proxy; secant update refines)

    if debug
        fprintf('switch to Broyden (bkw): normdr2=%g\n', norm(dr2_cur));
    end

    %% ====================================================================
    %% FORWARD BRANCH (mirror of prtlambertTB)
    %% ====================================================================
else
    v1_cur   = v1tnew;
    xout     = xoutfwd;
    dr2_cur  = xoutfwd(end,1:3) - r2;
    deltar2  = dr2_cur * Learning_rate;
    normdr2  = norm(deltar2);
    prop_r2  = r2;
    J          = eye(3);
    v1_t_prev  = v1tnew;
    dr2_t_prev = dr2_cur;

    n_thompson   = 0;
    normdr2_prev = normdr2;
    while normdr2 > T_switch && (n_thompson < N_switch || normdr2 > normdr2_prev) && integration_number < maxIter
        prop_r2  = prop_r2 - deltar2;
        [v1_cur, v2tnew_cur] = lambertb(r1, prop_r2, v1_cur, dm, de, nrev, dtsec,mu);
        if ~all(isreal(v1_cur)) || any(isnan(v1_cur))
            if debug,  history.converged = false;  end
             v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
            hitrad = NaN;  hitearthvar = NaN;
            return
        end
        if  hitearth ( 100, r1, v1_cur, r2, v2tnew_cur, nrev )
             v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
            errorb = 'Hit earth';
            return
        end
        oe = kep_elements(r1, v1_cur, mu);
        a  = oe(1);
        if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
            T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_orb = sqrt(a^3/mu) * 2*pi;
        end
        options_fwd.Sec = T_orb / options.delta;
        [~, xout]  = odeMPCI(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
        integration_number = integration_number + 1;
        if ~all(isfinite(xout(end,:))), break, end
        normdr2_prev = normdr2;
        deltar2   = (xout(end,1:3) - r2) * Learning_rate;
        normdr2   = norm(deltar2);
        dr2_t_cur = xout(end,1:3) - r2;
        dx_t = (v1_cur - v1_t_prev)';
        df_t = (dr2_t_cur - dr2_t_prev)';
        if norm(dx_t) > 0
            J = J + (df_t - J*dx_t) * dx_t' / (dx_t'*dx_t);
        end
        v1_t_prev  = v1_cur;
        dr2_t_prev = dr2_t_cur;
        n_thompson = n_thompson + 1;

        if debug
            history.iterations(end+1) = integration_number;
            history.normdr2(end+1)    = normdr2 / Learning_rate;
            history.phase{end+1}      = 'Thompson_fwd';
            fprintf('integration %d (T fwd): normdr2=%g\n', integration_number, normdr2/Learning_rate);
        end
    end
    v1_B  = v1_cur;
    dr2_B = xout(end,1:3) - r2;
    dr2_cur = dr2_B;

    if debug
        fprintf('switch to Broyden (fwd): normdr2=%g\n', norm(dr2_cur));
    end
end

%% Check if Thompson already converged
if norm(dr2_cur) <= Tolr
    if debug
        history.converged   = true;
        history.final_error = norm(dr2_cur);
    end
    v1t = v1_cur;
    v2t = xout(end,4:6);
    if  hitearth ( 100, r1, v1t, r2, v2t, nrev )
        v1t = zeros(1,3)+1i;  v2t = zeros(1,3)+1i;
        errorb = 'Hit earth';
    end
    return
end

%% ---- Broyden initialization (J from Phase-1; refine with lambertb→exit secant) ----
dx0 = (v1_B - v1_A)';
df0 = (dr2_B - dr2_A)';
if norm(dx0) > 0
    J = J + (df0 - J*dx0) * dx0' / (dx0'*dx0);
end

v1_cur = v1_B;

%% ---- Phase 2: Broyden forward ----
while norm(dr2_cur) > Tolr && integration_number < maxIter
    v1_next  = v1_cur - (J \ dr2_cur')';
    v1_prev  = v1_cur;
    dr2_prev = dr2_cur;
    v1_cur   = v1_next;

    oe = kep_elements(r1, v1_cur, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_fwd.Sec = T_orb / options.delta;

    [~, xout]  = odeMPCI(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    dr2_cur = xout(end,1:3) - r2;

    dx = (v1_cur  - v1_prev)';
    df = (dr2_cur - dr2_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if debug
        history.iterations(end+1) = integration_number;
        history.normdr2(end+1)    = norm(dr2_cur);
        history.phase{end+1}      = 'Broyden';
        fprintf('integration %d (Broyden): normdr2=%g\n', integration_number, norm(dr2_cur));
    end

    if ~all(isreal(v1_cur)) || any(~isfinite(v1_cur)) || ~isfinite(norm(dr2_cur))
        break
    end
end

%% Convergence check
if ~isfinite(norm(dr2_cur)) || norm(dr2_cur) > Tolr || hitearth(100, r1, v1_cur, r2, xout(end,4:6), nrev)
    %warning('prtlambertTRB: did not converge');
    if debug
        history.converged   = false;
        history.final_error = NaN;
    end
    v1t = v1_cur + 1i;
    v2t = xout(end,4:6) + 1i;
else
    v1t = v1_cur;
    v2t = xout(end,4:6);
end

if debug
    history.converged   = true;
    history.final_error = norm(dr2_cur);
end



%[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
end
