function [Y_pop,IDX] = EvaluateModel(X,Params,model)
%EVALUATEMODEL  Per-leg transfer cost function for the lower level.
%   [Y_pop,IDX] = EVALUATEMODEL(X, Params, model) computes the total DeltaV of each
%   candidate transfer in X by solving the perturbed Lambert problem (prtlambertT)
%   against the target propagated by its Chebyshev model. Used by the one-variable
%   lower-level optimizer.
%
%   Inputs:
%     X      - candidate matrix: column 1 = normalized transfer time,
%              column 2 = normalized discrete-setting index
%     Params - struct of packed leg state + environment (rc, vc, T, rt1, vt1, T1,
%              mc, Ac, J, Cd, Re, mu1, RelativeEpoch)
%     model  - target Chebyshev propagation model
%   Outputs:
%     Y_pop  - total DeltaV per candidate [km/s] (1000 flags infeasible/failed)
%     IDX    - chosen discrete setting index (1..16) per candidate

T = Params.T;
RelativeEpoch = Params.RelativeEpoch;
t1 = Params.T1;
rt1 = Params.rt1;
vt1 = Params.vt1;
rc = Params.rc;
vc = Params.vc;
mu1 = Params.mu1;
Cd = Params.Cd;
Ac = Params.Ac;
mc = Params.mc;
Re = Params.Re;
J = Params.J;
[at,~,~,~,~,~] = kep_elements(rt1,vt1,mu1);
Torbt = sqrt(at^3/mu1)*2*pi;
[ac,~,~,~,~,~] = kep_elements(rc,vc,mu1);
Torbc = sqrt(ac^3/mu1)*2*pi;
options.AbsTol = 1e-7;
options.RelTol = 1e-5;
options.N = 16;
delta = 8;
options.Sec = min(Torbt,Torbc)/delta;

Y_pop = nan(size(X,1),1);
IDX = nan(size(X,1),1);
TM = X(:,1);
Sett = round(X(:,2)*15)+1;
%NREV = X(:,1);
parfor N = 1 : size(X,1)
    %nrev = round(NREV(N)*3,0);
    tm = TM(N)*3*60*60;
    I = Sett(N);
    if RelativeEpoch > 2*60*60*24*30*12 - tm - t1 - T
        Y_pop(N,1) = 1000;
    else
        xt =  propCheb(RelativeEpoch+1e-10 + t1 + T + tm,model);
        rt2 = xt(1:3);
        vt2 = xt(4:6);
        %Y = inf*ones(16,1);
        DM = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
        DE = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
        Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];
        
            [vc1m,vc2m] = prtlambertT(@(t,x)orbit_eq_J6_drag(t,x,mu1,Cd,Ac,mc,Re,J),@odeMPCI,rc,rt2,vc,DM(I),DE(I),Nrev(I), tm ,1,options);
            if isreal(vc1m) && isreal(vc2m)
                deltav1 = norm(vc1m-vc);
                deltav2 = norm(vc2m-vt2);
                deltavtot = deltav1+deltav2;

            else
                deltavtot = 1000;
            end
            %Y(I) = deltavtot;
        
        %[M, Idx] = min(Y);
        Y_pop(N,1) = deltavtot;
        IDX(N,1) = I;
    end
end

end

