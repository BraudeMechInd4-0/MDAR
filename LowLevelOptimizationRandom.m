function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationRandom(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1,Solver)
%LOWLEVELOPTIMIZATIONRANDOM  Lower-level transfer optimizer (random search).
%   [DeltaV,tm,de,dm,nrev] = LOWLEVELOPTIMIZATIONRANDOM(rc, vc, T, RelativeEpoch,
%   model, rt1, vt1, T1, mc, Ac, J, Cd, Re, mu1, Solver) samples random transfer
%   parameters (normalized time + a 4-bit direction/energy/revolution setting),
%   evaluates each with the perturbed-Lambert lower level, and returns the
%   minimum-DeltaV maneuver.
%
%   Inputs:
%     rc, vc        - current chaser position/velocity
%     T             - elapsed mission time [s]
%     RelativeEpoch - target's relative epoch offset [s]
%     model         - target Chebyshev propagation model
%     rt1, vt1      - target state at reference epoch
%     T1            - planned wait before the maneuver [s]
%     mc, Ac, J, Cd, Re, mu1 - chaser/environment parameters
%     Solver        - Lambert-solver select (with the J6+Drag+SRP+Moon force model):
%                     1 -> prtlambertT (default), 2 -> prtlambertTRB,
%                     3 -> prtlambertMPS, 4 -> prtlambertMPSRB
%   Outputs:
%     DeltaV        - best total DeltaV found among samples [km/s]
%     tm            - corresponding transfer time [s]
%     de, dm        - discrete energy / direction-of-motion flags
%     nrev          - number of revolutions
nVar = 2;                % Number of Unknown (Decision) Variables
%VarSize=[1 nVar];       
VarMin = 1e-5;             
VarMax = 1;             
Params.rc = rc;
Params.vc = vc;
Params.T = T;
Params.RelativeEpoch = RelativeEpoch;
Params.rt1 = rt1;
Params.vt1 = vt1;
Params.T1 = T1;
Params.Ac = Ac;
Params.mc = mc;
Params.J = J;
Params.mu1 = mu1;
Params.Cd = Cd;
Params.Re = Re;

% Maximum Number of Iterations
MaxIt = 20;

% Population Size (and Number of Offsprings)
lambda = (4+round(3*log(nVar)))*2;



%Overall number of solutions 
Nsolutions = MaxIt*(lambda + 1) + 1;

% Generate solutions

X = [rand(Nsolutions,1)];
Idx = randi([0,1],Nsolutions,4);



[Cost,IDX]=EvaluateModel(X,Idx,Params,model,Solver);

[DeltaV, I] = min(Cost);


DM = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
DE = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];
tm = X(I)*4*60*60;
de = DE(IDX(I));
dm = DM(IDX(I));
nrev = Nrev(IDX(I));

end

function [Y_pop,IDX] = EvaluateModel(X,Idx,Params,model,Solver)

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
Sett = sum(Idx.*[1 2 4 8],2)+1;
%NREV = X(:,1);
parfor N = 1 : size(X,1)
    %nrev = round(NREV(N)*3,0);
    tm = TM(N)*4*60*60;
    I = Sett(N);
    if RelativeEpoch > 2*60*60*24*30*12 - tm - t1 - T
        Y_pop(N,1) = 1000;
        IDX(N,1) = I;
    else
        xt =  propCheb(RelativeEpoch+1e-10 + t1 + T + tm,model);
        rt2 = xt(1:3);
        vt2 = xt(4:6);
        %Y = inf*ones(16,1);
        DM = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
        DE = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
        Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];
        if Solver == 2
            [vc1m,vc2m] = prtlambertTRB(@(t,x)orbit_eq_J6_drag_SRP_moon(t,x,mu1,Cd,Ac,mc,Re,J),rc,rt2,vc,DM(I),DE(I),Nrev(I), tm ,1,options);
                            
        elseif Solver == 3
            [vc1m,vc2m,~] = prtlambertMPS(@(t,x)orbit_eq_J6_drag_SRP_moon(t,x,mu1,Cd,Ac,mc,Re,J),@odeMPCI,rc,rt2,vc,DM(I),DE(I),Nrev(I), tm ,1e-5,options);
        elseif Solver == 4
            [vc1m,vc2m,~] = prtlambertMPSRB(@(t,x)orbit_eq_J6_drag_SRP_moon(t,x,mu1,Cd,Ac,mc,Re,J),@odeMPCI,rc,rt2,vc,DM(I),DE(I),Nrev(I), tm ,1e-5,options);
        else
            [vc1m,vc2m] = prtlambertT(@(t,x)orbit_eq_J6_drag_SRP_moon(t,x,mu1,Cd,Ac,mc,Re,J),@odeMPCI,rc,rt2,vc,DM(I),DE(I),Nrev(I), tm ,1,options);
                                        
        end
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

