function [tout,xout,errorOutput,model] = odeMPCI(ODEFUN, TSPAN, X0,options)
%[tout,xout,error] = odeMPCI(ODEFUN, TSPAN,X0)
%[tout,xout,error] = odeMPCI(ODEFUN, TSPAN,X0,options)
%[tout,xout,error,model] = odeMPCI(ODEFUN, TSPAN, X0,options)
%Bai's [1] MPCI integration method, as described in [2]
% Useage
%       ODEFUN - the function to integrate (the f in the equation dx/dt = f(x,t))
%       TSPAN  - The time span of the integration. Can be [Tstart,Tend] or a
%                vector of timepoints to be evaluated
%       X0     - Initial conditions (shoulve be given as x(t = 0))
%       options.AbsTol - Absolute tolerance
%       options.RelTol - Relative tolerace
%       options.N      - number of timepoints to be sampled. Optional if TSPAN is a
%                 vector longer than 2, mandatory if it is [Tstart,Tend]
%       options.deltaT - The timestep over each to fit the polynomial.
%       options.xtinit - An initial guess. For astrodynamics it is customery to use
%                the solution for the two body problem, or even better one
%                of the SGP analytical methods.
%       tout   - A list of time points
%       xout   - Values of the function approximation on these points
%       error  - error code: 0 if all went well. -1 if failed to converge
%       model - a structure containing the chebychev polynomials.
%       model.ts - a list of timestamps, where each polynomial starts
%       model.b  - a matrix, each line is a set of polynomial coefficients
%                  the first line matches the first time in ts. The first
%                  coefficient is b_0
%                  To compute the function at time t you need sum(T.*b)
%[1] Bai, X. and Junkins, J.L., 2011. Modified Chebyshev-Picard iteration methods for solution of boundary value problems. The Journal of the astronautical sciences, 58(4), pp.615-642.
%[2] Woollands, R. and Junkins, J.L., 2019. Nonlinear differential equation solvers via adaptive Picard–Chebyshev iteration: Applications in astrodynamics. Journal of Guidance, Control, and Dynamics, 42(5), pp.1007-1022.
%See also:  Darin Koblick (2024). Vectorized Picard-Chebyshev Method (https://www.mathworks.com/matlabcentral/fileexchange/36940-vectorized-picard-chebyshev-method), MATLAB Central File Exchange. Retrieved October, 2023. 

%% input check and initialization
errorOutput = 0;
if nargin <3
    error("not enough input arguments");
end
if ~isfield(options, 'AbsTol')
    AbsTol = 1e-12;
else
    AbsTol = options.AbsTol;
end
if ~isfield(options, 'RelTol')
    RelTol = 1e-9;
else
    RelTol = options.AbsTol;
end
if  ~isfield(options, 'N')
    if length(TSPAN) > 2
        N = length(TSPAN);
    else
        N = 32;
    end
else
    N = options.N;
end
if ~isfield(options,'Sec') || isnan(options.Sec) || ~isreal(options.Sec)
    Sec = TSPAN(end);
else
    Sec= options.Sec;
end
outputmodel = false;
if nargout ==4
    outputmodel = true;
end
N = N-1; %because it's 0 - N in all places so if N is 32 and we don't deduct one we get 33 points. 
tstart = 0;
tend = Sec;
output_length = N*ceil(TSPAN(end)/Sec)+1;
firstpos = 1;
lastpos = NaN;
xout = NaN(output_length,6);
xout(1,:) = X0;
tout = NaN(output_length,1);

if outputmodel
    model_length = floor(TSPAN(end)/Sec)+1;
    model_counter = 1;
    model.ts = zeros(model_length+1,1);
    model.ts(1) = 0;
    model.bi = NaN(model_length,N+1,6);
else
    model = [];
end

tau = fliplr(cos((0:N)*(pi/N))); %Gauss lobato nodes
if ~isfield(options,'xinit')
    [m,n] = size(X0);
    if m > 1 && n >=1 %this means it's a matrix, so we have different X's so X0 is not the boundery value but an initial guess.
        error("if you want an initial guess use xtinit")
    end
    if n == 1 && m ~=1 %it is a column instead of a vector.
        X0 = X0';
    end

    xtinit = repmat(X0,[length(tau),1]); %innitial guess
end

X0 = [X0;zeros(N,length(X0))];
lastrunflag = false;
while true
    xold = xtinit;%zeros(length(tau),length(X0)); %innitial guess
    om2 = (tend-tstart)/2;
    om1 = (tend+tstart)/2;

    %% First Creat the Vectors and matrices
    W = eye(length(tau));
    W(1,1) = 0.5;
    W(end,end) = 0.5;
    T = cos((0:N-1)'*acos(tau))'; %eq A6 k = 0,1,2,...N-1
    Tm1 = cos((0:N).*pi()); %acos(-1) = pi;
    L = [Tm1;zeros(N,N+1)];
    s_ = 1./(4:2:2*N);
    S_3 = zeros(N+1);
    S_3(2:end,2:end) = diag([-.5,-s_(1:end-1)],1);
    S_2 = diag([1,s_],-1);
    S_1 = S_2+S_3;
    S = S_1(:,1:N);
    S(1,:) = [1/4, zeros(1,N-1)];
    A = (T'*W*T)\T'*W;


    %% Get the function g using picard iterations

    T = cos((0:N)'*acos(tau))';   
    eAbs = inf;
    eRel = inf;
    i = 0;
    while (eAbs  > AbsTol || eRel > RelTol) && i<2000
        F = ODEFUN(om2.*tau+om1,xold); 
        P1 = om2*(eye(N+1) - L) * S;
        bi = X0 + P1 * A * F;
        xnew = T*bi;
        eAbs = max(abs(xnew - xold),[],'all');
        eRel = max(abs(xnew - xold)./min(xnew,xold),[],'all');
        xold = xnew;
        i = i+1;
    end

    if i >= 4000
        %error("did not converge")
        errorOutput = -1;
        tout = 0;
        xout=0;
        return
    end


    %% now evaluate at the required t's
    %for each t do x(t) = sum(beta_k*T_k(t))
    
    if length(TSPAN) ~= 2
        %evaluate at the points in TSPAN and return
        tau = -(TSPAN(end)-TSPAN(1))/(TSPAN(end)+TSPAN(1)) +2*(TSPAN)/(TSPAN(end)-TSPAN(1));
        tout = TSPAN;
        T = cos((0:N)'.*acos(tau))'; %length(tau)xN+1
        xout_ = T * bi;
    else
        lastpos = firstpos+N;
        tout(firstpos:lastpos) = tstart + (tau+1)*(tend-tstart)/2;
        xout_ = xnew;
    end
    
    
    xout(firstpos+1:lastpos,:) = xout_(2:end,:);
    if outputmodel
        model.bi(model_counter,:,:) = bi;
        model_counter = model_counter+1;
        model.ts(model_counter) = tend;
    end
    if lastrunflag
        return
    end
    firstpos = lastpos;
    tstart = tend;
    tend = tstart+Sec;
    
    if tend > TSPAN(end)
        tend = TSPAN(end);
        lastrunflag = true;
    end
    X0 = xout_(end,:);
    xtinit = repmat(X0,[length(tau),1]); %innitial guess
    X0 = [X0;zeros(N,length(X0))];
end
