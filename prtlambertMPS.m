function [v1t, v2t, hitearthvar, hitrad, integration_number, history] = prtlambertMPS(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTMPS  Perturbed Lambert solver — pure Method of Particular Solutions (forward only).
%
% Usage:
%   [v1t,v2t,hitearthvar,hitrad,integration_number,history] = prtlambertMPS( ...
%       forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Inputs:
%   forcemodel   - function handle @(t,x) for the perturbing force model
%   odesolver    - function handle to the integrator, e.g. @odeMPCI
%   options      - struct with fields: maxIter, Tolr, delta, N, debug
%   r1           - initial position (1x3, km)
%   r2           - target  position (1x3, km)
%   v1           - reference velocity for lambertb seed (1x3, km/s)
%   dm           - direction of motion ('L' or 'S')
%   de           - orbital energy    ('L' or 'H')
%   nrev         - number of revolutions
%   dtsec        - transfer time (s)
%
% Outputs:
%   v1t              - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t              - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   hitearthvar      - hitearth flag
%   hitrad           - closest approach radius (km)
%   integration_number - total number of integrations performed
%   history          - struct with convergence data (populated only when options.debug=true)
%                      Fields: iterations, integration_number, normdr2, converged, final_error

mu = 398600.5;
r_GEO = 42164;
r_max_allowed = 1.5 * r_GEO;

if nargin < 10,  options = struct();  end
if ~isfield(options, 'Tolr'),    options.Tolr    = 1e-8;  end
if ~isfield(options, 'maxIter'), options.maxIter = 100;   end
if ~isfield(options, 'delta'),   options.delta   = 16;    end

debug = isfield(options, 'debug') && options.debug;
if debug
    fprintf('starting MPS\n');
    history = struct('iterations', [], 'integration_number', [], 'normdr2', [], ...
                     'converged', false, 'final_error', NaN);
else
    history = struct();
end

[v1ref] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1ref)) || any(isnan(v1ref)) || any(~isfinite(v1ref))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    hitrad = NaN;  hitearthvar = NaN;
    integration_number = 0;
    return
end

xref   = [r1 v1ref];
oeref  = kep_elements(r1, v1ref, mu);
a_ref  = oeref(1);
if a_ref <= 0 || ~isreal(a_ref) || isnan(a_ref) || a_ref > 2*r_max_allowed
    Torbref = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    Torbref = sqrt(a_ref^3/mu) * 2*pi;
end
options.Sec = Torbref / options.delta;

delta_size = 1e-7;
i = 1;
integration_number = 0;
Tolr    = options.Tolr;
maxIter = floor(options.maxIter / 4);   % 4 integrations per MPS iteration

%% Pre-loop: integrate nominal, compute initial rdep
[~, xnew] = odesolver(forcemodel, [0, dtsec], xref, options);
integration_number = integration_number + 1;
rdep = r2 - xnew(end, 1:3);

if debug
    history.iterations(end+1)         = 0;
    history.integration_number(end+1) = integration_number;
    history.normdr2(end+1)            = norm(rdep);
end

%% Main loop
xpert = zeros(3,6);
while norm(rdep) > Tolr && i < maxIter
    normv0   = norm(xref(4:6));
    deltardot = [zeros(3), eye(3)*delta_size*normv0];
    xprtrbd  = deltardot + xref;

    for j = 1:3
        [~, xtmp]   = odesolver(forcemodel, [0, dtsec], xprtrbd(j,:), options);
        xpert(j,:)  = xtmp(end,:);
    end
    integration_number = integration_number + 3;
    if ~all(isfinite(xpert(:))), break, end

    if debug
        fprintf('iteration %d: error = %g km\n', i, norm(rdep));
    end

    DR = (xpert(:,1:3) - xnew(end,1:3))';
    if rcond(DR) < 1e-10
        warning('prtlambertMPS: DR near-singular (rcond=%.2e), stopping', rcond(DR));
        break;
    end
    alpha = (DR \ rdep')';

    xref  = xref + alpha * deltardot;
    oeref = kep_elements(xref(1:3), xref(4:6), mu);
    a_ref = oeref(1);
    if a_ref <= 0 || ~isreal(a_ref) || isnan(a_ref) || a_ref > 2*r_max_allowed
        Torbref = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        Torbref = sqrt(a_ref^3/mu) * 2*pi;
    end
    options.Sec = Torbref / options.delta;

    [~, xnew] = odesolver(forcemodel, [0, dtsec], xref, options);
    integration_number = integration_number + 1;
    if ~all(isfinite(xnew(end,:))), break, end
    rdep = r2 - xnew(end, 1:3);

    if debug
        history.iterations(end+1)         = i;
        history.integration_number(end+1) = integration_number;
        history.normdr2(end+1)            = norm(rdep);
    end
    i = i + 1;
end

if i >= maxIter || ~isfinite(norm(rdep)) || norm(rdep) > Tolr
    warning('prtlambertMPS: did not converge');
    if debug
        history.converged   = false;
        history.final_error = NaN;
    end
    v1t = NaN(1,3);  v2t = NaN(1,3);
    hitearthvar = NaN;  hitrad = NaN;
else
    if debug
        history.converged   = true;
        history.final_error = norm(rdep);
    end
    v1t = xnew(1,   4:6);
    v2t = xnew(end, 4:6);
    [hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
end
end
