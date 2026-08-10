function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationGA(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1)
%LOWLEVELOPTIMIZATIONGA  Lower-level transfer optimizer (genetic algorithm).
%   [DeltaV,tm,de,dm,nrev] = LOWLEVELOPTIMIZATIONGA(rc, vc, T, RelativeEpoch, model,
%   rt1, vt1, T1, mc, Ac, J, Cd, Re, mu1) searches a single transfer's parameters
%   with a GA and returns the best maneuver found via the perturbed-Lambert solver.
%
%   Inputs:
%     rc, vc        - current chaser position/velocity
%     T             - elapsed mission time [s]
%     RelativeEpoch - target's relative epoch offset [s]
%     model         - target Chebyshev propagation model
%     rt1, vt1      - target state at reference epoch
%     T1            - planned wait before the maneuver [s]
%     mc, Ac, J, Cd, Re, mu1 - chaser/environment parameters
%   Outputs:
%     DeltaV        - best total DeltaV found [km/s]
%     tm            - transfer time [s]
%     de, dm        - discrete energy / direction-of-motion flags
%     nrev          - number of revolutions
nVar = 2;                % Number of Unknown (Decision) Variables
%VarSize=[1 nVar];
VarMin = 10e-5;
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
PopSize = (4+round(3*log(nVar)))*2 + 1;


% Generate solutions
ParentsX = [rand(PopSize,1)];
ParentsI = randi([0,1],PopSize,4);


[CostParent,IDX]=EvaluateModel(ParentsX,ParentsI,Params,model);

[DeltaV, I] = min(CostParent);
Xbest = ParentsX(I);
IDXBest = IDX(I);

for gen = 1 : MaxIt
    [MTX,MTI] = SelectionLow(ParentsX,ParentsI,CostParent,PopSize);
    [OffspringX,OffspringI] = ReproductionLow(MTX,MTI,PopSize);
    [CostOffspring,IDX]=EvaluateModel(OffspringX,OffspringI,Params,model);
    if min(CostOffspring) < DeltaV
        [DeltaV, I] = min(CostOffspring);
        Xbest = OffspringX(I);
        IDXBest = IDX(I);
    end

    [ParentsX,ParentsI,CostParent] = updatePop(ParentsX,ParentsI,CostParent,OffspringX,OffspringI,CostOffspring,PopSize);
end
DM = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
DE = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];
tm = Xbest*3*60*60;
de = DE(IDXBest);
dm = DM(IDXBest);
nrev = Nrev(IDXBest);

end

function [Y_pop,IDX] = EvaluateModel(X,Idx,Params,model)

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


function [MTX,MTI] = SelectionLow(ParentsX,ParentsI,CostParent,PopSize)
MTX = nan(PopSize,1);
MTI = nan(PopSize,4);
for i = 1 : PopSize
    S = randi(PopSize,2,1);
    if CostParent(S(1)) < CostParent(S(1))
        MTX(i,:) = ParentsX(S(1),:);
        MTI(i,:) = ParentsI(S(1),:);
    else
        MTX(i,:) = ParentsX(S(2),:);
        MTI(i,:) = ParentsI(S(2),:);
    end
end
end

function [OffspringX,OffspringI] = ReproductionLow(MTX,MTI,PopSize)
Offspring1X = MTX;
Offspring1I = MTI;
%% Crossover
pc = 0.8;
for i = 1 : 2 :PopSize-1
    if rand() < pc
        alpha = rand();
        Offspring1X(i,:) = MTX(i,:)*alpha + MTX(i+1,:)*(1-alpha);
        Offspring1X(i+1,:) = MTX(i+1,:)*alpha + MTX(i,:)*(1-alpha);
        point = randi(length(Offspring1I(i,:))-1);
        Offspring1I(i,:) =  [MTI(i,1:point),MTI(i+1,point+1:end)];
        Offspring1I(i+1,:) =  [MTI(i+1,1:point),MTI(i,point+1:end)];
    end
end
%Mutation
pm = 0.1;
OffspringX = Offspring1X;
OffspringI = Offspring1I;
for i = 1 : PopSize
    if rand < pm
        OffspringX(i,:) = rand(size(Offspring1X(i,:)));
        point = randi(length(OffspringI(i,:)));
        OffspringI(i,point) = 1 - OffspringI(i,point);
    end
end
end

function [NextParentsX,NextParentsI,NextParentCost] = updatePop(ParentsX,ParentsI,CostParent,OffspringX,OffspringI,CostOffspring,PopSize)
PopAllX = [ParentsX;OffspringX];
PopAllI = [ParentsI;OffspringI];
CostAll = [CostParent;CostOffspring];
[~, idx] = sort(CostAll);
NextParentsX = PopAllX(idx(1:PopSize),:);
NextParentsI = PopAllI(idx(1:PopSize),:);
NextParentCost = CostAll(idx(1:PopSize));
end
