function [v1t, v2t, errorb] = prtlambertbf(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, mu1, Tolr, solveroptions)
%
% Solve the perturbed Lambert's problem using Thompson iteration.
%Usage:
%    [v1t, v2t, errorb] = prtlambertbf(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, mu1, Tolr, solveroptions)
%
% Inputs:
%   forcemodel   - function handle: @(t,x) orbit_eq_J6_drag(t,x,mu,Cd,A,m,Re,J)
%   odesolver    - function handle: @odeMPCI
%   r1           - initial position [km]
%   r2           - target  position [km]
%   v1           - initial velocity [km/s]
%   dm           - direction of motion ('L' or 'S')
%   de           - orbital energy    ('L' or 'H')
%   nrev         - number of revolutions
%   dtsec        - transfer time [s]
%   mu1          - gravitational parameter [km^3/s^2]
%   Tolr         - position convergence tolerance [km]
%   solveroptions - struct passed to odesolver; may contain:
%                   .maxIter       (default 20)
%                   .Learning_rate (default 1)
%                   .delta         (default 8, ODE step = T_orb/delta)
%                   .debug         (default false)
%
% Outputs:
%   v1t    - transfer velocity at r1 [km/s]  (complex if solver failed)
%   v2t    - transfer velocity at r2 [km/s]  (complex if solver failed)
%   errorb - hitearth flag (from hitearth()); NaN if solver failed

r_GEO         = 42164;           % GEO radius [km]
r_max_allowed = 1.5 * r_GEO;    % cap for hyperbolic/very-large orbits

% Read optional fields from solveroptions
if isfield(solveroptions, 'maxIter')
    maxIter = solveroptions.maxIter;
else
    maxIter = 20;
end
if isfield(solveroptions, 'Learning_rate')
    Learning_rate = solveroptions.Learning_rate;
else
    Learning_rate = 1;
end
if isfield(solveroptions, 'delta')
    delta = solveroptions.delta;
else
    delta = 8;
end
debug = isfield(solveroptions, 'debug') && solveroptions.debug;

prop_r2 = r2;
deltar  = inf;
i       = 1;

while norm(deltar) >= Tolr && i < maxIter

    [v1t, v2t] = lambertb(r1, prop_r2, v1, dm, de, nrev, dtsec);

    if ~all(isreal(v1t)) || ~all(isreal(v2t)) || all(isnan(v1t)) || all(isnan(v2t))
        if debug
            fprintf("prtlambertbf: Lambert failed at iteration %d\n", i);
        end
        v1t    = v1t + 1i;
        v2t    = v2t + 1i;
        errorb = NaN;
        return
    end

    % Compute adaptive ODE step from current orbit period
    oe = kep_elements(r1, v1t, mu1);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi * sqrt((2*r_max_allowed)^3 / mu1);
    else
        T_orb = 2*pi * sqrt(a^3 / mu1);
    end
    solveroptions.Sec = T_orb / delta;

    [~, x2] = odesolver(forcemodel, [0 dtsec], [r1 v1t], solveroptions);
    r_new = x2(end, 1:3);
    v2t   = x2(end, 4:6);

    deltar  = (r_new - r2) * Learning_rate;
    prop_r2 = prop_r2 - deltar;

    if debug
        fprintf("prtlambertbf: iteration %d, error = %g km\n", i, norm(deltar));
    end

    i = i + 1;
end

if i >= maxIter || ~isfinite(norm(deltar)) || norm(deltar) > Tolr
    warning("prtlambertbf: Thompson did not converge after %d iterations", maxIter);
    v1t    = v1t + 1i;
    v2t    = v2t + 1i;
    errorb = NaN;
    return
end

[errorb, ~] = hitearth(100, r1, v1t, r2, v2t, nrev);
