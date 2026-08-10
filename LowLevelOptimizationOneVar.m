function [DeltaV,tm,de,dm,nrev] = LowLevelOptimizationOneVar(rc,vc,T,RelativeEpoch,model,rt1,vt1,T1,mc,Ac,J,Cd,Re,mu1)
%LOWLEVELOPTIMIZATIONONEVAR  Lower-level transfer optimizer (CMA-ES).
%   [DeltaV,tm,de,dm,nrev] = LOWLEVELOPTIMIZATIONONEVAR(rc, vc, T, RelativeEpoch,
%   model, rt1, vt1, T1, mc, Ac, J, Cd, Re, mu1) optimizes a single transfer with a
%   CMA-ES search over the normalized transfer time and discrete mode, scoring
%   candidates via the perturbed-Lambert lower level (EvaluateModel).
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

VarSize=[1 nVar];       % Decision Variables Matrix Size

VarMin=0;             % Lower Bound of Decision Variables
VarMax= 1;             % Upper Bound of Decision Variables
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

%% CMA-ES Settings

% Maximum Number of Iterations
MaxIt = 10;

% Population Size (and Number of Offsprings)
lambda = 2*(4+round(3*log(nVar)))*2;

% Number of Parents
mu=round(lambda/2);

% Parent Weights
w=log(mu+0.5)-log(1:mu);
w=w/sum(w);

% Number of Effective Solutions
mu_eff=1/sum(w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
sigma0=0.3*(VarMax-VarMin);
cs=(mu_eff+2)/(nVar+mu_eff+5);
ds=1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Covariance Update Parameters
cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1=2/((nVar+1.3)^2+mu_eff);
alpha_mu=2;
cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth=(1.4+2/(nVar+1))*ENN;


ps=cell(MaxIt,1);
pc=cell(MaxIt,1);
C=cell(MaxIt,1);
sigma=cell(MaxIt,1);

ps{1}=zeros(VarSize);
pc{1}=zeros(VarSize);
C{1}=eye(nVar);
sigma{1}=sigma0;

empty_individual.Position=[];
empty_individual.Step=[];
empty_individual.Cost=[];

M=repmat(empty_individual,MaxIt,1);
M(1).Position=unifrnd(VarMin,VarMax,VarSize);
M(1).Step=zeros(VarSize);
[M(1).Cost,M(1).IDX]=EvaluateModel(M(1).Position,Params,model);

BestSol=M(1);

BestCost=zeros(MaxIt,1);

%% CMA-ES Main Loop

for g=1:MaxIt

    % Generate Samples
    pop=repmat(empty_individual,lambda,1);
    X = nan(lambda,nVar);
    for i=1:lambda
        pop(i).Step=mvnrnd(zeros(VarSize),C{g});
        pop(i).Position = M(g).Position+sigma{g}*pop(i).Step;
        pop(i).Position = max([pop(i).Position;1e-10*ones(size(pop(i).Position))],[],1);
        pop(i).Position = min([pop(i).Position;ones(size(pop(i).Position))],[],1);
        X(i,:) = pop(i).Position;
    end

    [CurrentCosts,CurrentIDX]=EvaluateModel(X,Params,model);
    for i=1:lambda
        % Update Best Solution Ever Found
        pop(i).Cost = CurrentCosts(i);
        pop(i).IDX = CurrentIDX(i);
        if pop(i).Cost<BestSol.Cost
            BestSol=pop(i);
        end
    end
    % Sort Population
    Costs=[pop.Cost];
    [~, SortOrder]=sort(Costs);
    pop=pop(SortOrder);

    % Save Results
    BestCost(g)=BestSol.Cost;

    % Display Results
    %disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);

    % Exit At Last Iteration
    if g==MaxIt
        break;
    end

    % Update Mean
    M(g+1).Step=0;
    for j=1:mu
        M(g+1).Step=M(g+1).Step+w(j)*pop(j).Step;
    end
    M(g+1).Position=M(g).Position+sigma{g}*M(g+1).Step;
    M(g+1).Position = max([M(g+1).Position;1e-20*ones(size(M(g+1).Position))],[],1);
    M(g+1).Position = min([M(g+1).Position;ones(size(M(g+1).Position))],[],1);
    [M(g+1).Cost,M(g+1).IDX]=EvaluateModel(M(g+1).Position,Params,model);
    if M(g+1).Cost<BestSol.Cost
        BestSol=M(g+1);
    end

    % Update Step Size
    ps{g+1}=(1-cs)*ps{g}+sqrt(cs*(2-cs)*mu_eff)*M(g+1).Step/chol(C{g})';
    sigma{g+1}=sigma{g}*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;

    % Update Covariance Matrix
    if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1)))<hth
        hs=1;
    else
        hs=0;
    end
    delta=(1-hs)*cc*(2-cc);
    pc{g+1}=(1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
    C{g+1}=(1-c1-cmu)*C{g}+c1*(pc{g+1}'*pc{g+1}+delta*C{g});
    for j=1:mu
        C{g+1}=C{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
    end

    % If Covariance Matrix is not Positive Defenite or Near Singular
    [V, E]=eig(C{g+1});
    if any(diag(E)<0)
        E=max(E,0);
        C{g+1}=V*E/V;
    end

end
DeltaV = BestSol.Cost;
DM = ['L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S'];
DE = ['H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L';'H';'H';'L';'L'];
Nrev = [0;0;0;0;0;1;1;1;1;2;2;2;2;3;3;3;3];
tm = BestSol.Position(1)*3*60*60;
de = DE(BestSol.IDX);
dm = DM(BestSol.IDX);
nrev = Nrev(BestSol.IDX);

end

