%HILEVEL  Single-objective bi-level MDAR optimization (paper 1).
%   Upper-level genetic algorithm that designs multi-debris removal missions,
%   evaluating each plan with a selectable lower-level trajectory optimizer.
%
%   Low-level type (set LowLevelType): "Rand" (random search), "GA", or "CMA".
%   Workflow:
%     1) Load Chebyshev models, constants, and the Iridium-33 debris catalog.
%     2) Initialize the upper-level population and evaluate it (EvaluateHi).
%     3) Iterate selection (tournamentSelection), crossover (CrossOver),
%        mutation (Mutation), evaluation, and elitism (EliteProcedure).
%     4) Save the best solution per run and low-level type.
%
%   Requires submodules/folders: ParseGP, Satellites, fundamentals-of-astrodynamics,
%   models. Loads: models/Iridium33Model.mat and Iridium33.mat.
%   Key parameters (below): PopSizeHi, MaxGenHi, NtoRemove, T1Max, weights wM/wDV.
clc;clear;close all
rng('shuffle')
%% Single-objective bi-level optimization run
% General Parameters + propagating the Satellite locaions

fprintf('loading...');
addpath(genpath("ParseGP"))
addpath(genpath("Satellites"))
addpath(genpath("fundamentals-of-astrodynamics"))
addpath(genpath("models"))
load("models\Iridium33Model.mat")
% if exist('mu')
%     path = which('mu');
%     rmpath(path(1:end-4));
% end

%%
% fprintf('loading...');
% if ~exist('consts','var')
%     [xpdotp, Re, J2, mu, whichconsts, consts, ~, dAT] = generate_parameters;
% end

%tumin = consts(1);
Cd = 2.2;
mu1 = 398600.5;
Re = 6378.137;
J = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, -0.2277161016e-6 ,0.5396484906e-6];
%xke    = 60.0 / sqrt(Re*Re*Re/mu);
%tumin  = 1.0 / xke;
%delta = 16;
%% generate list of objects and parameters: now it is for the chaser and one target
load('Iridium33.mat');
%save Satellites\IriduimList.mat m A epoch r v
fprintf("done!\n")

%RunIDX = 2;
PopSizeHi = 30;
MaxGenHi = 50;
mc = 500;%m(10);
Ac = 9e-6;%A(10);
NtoRemove = 10;
T1Max =  2*60*60*24*30;
GenParam.pm = 0.4;
GenParam.pc = 0.8;
wM = 1;
wDV = 8;

GenParam.PopSize = PopSizeHi;
GenParam.T1max = T1Max;
GenParam.NtoRemove = NtoRemove;
GenParam.NDeb = length(m);
GenParam.EliteMode = "H";
[~,ILastUpdated] = min(days(datetime("now") - epoch));
RelativeEpoch = seconds(epoch(ILastUpdated) - epoch);
Runs = 3;

mode = 1;Lstart = 2;
for RunIDX = 3:Runs
    for Ltype = 2 : 2
        
        if Ltype == 1
            LowLevelType = "Rand"; %Rand,GA,CMA
        elseif Ltype == 2
            LowLevelType = "GA"; %Rand,GA,CMA
        else
            LowLevelType = "CMA" ; %Rand,GA,CMA
        end
        fprintf('\nRun %d LowLeverType %s Start ... ', RunIDX, LowLevelType)
        Population = initPopHi(PopSizeHi,length(m),mode,NtoRemove,T1Max);
        %tic
        Population = EvaluateHi(Population,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType);
        %toc
        gen = 0;
        %BestSol = struct;
        [~,idxBest] = max([Population(:).fit]);
        BestSol = Population(idxBest);
        Ifmax = 1; % 1 for max, 0 for min
        %BestSol(MaxGenHi+1) = BestSol(1);
        fprintf('\n Generation %d Best perfprmance %4.3f with DV %f',gen, BestSol(1).Mremoved,sum(BestSol(1).DV))
        for gen = 1 : MaxGenHi
            PopulationSel = tournamentSelection(Population,Ifmax);
            PopulationCO = CrossOver(PopulationSel,GenParam);
            OffspringPop = Mutation(PopulationCO,GenParam);
            OffspringPop = EvaluateHi(OffspringPop,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType);
            [Population, BestSol(gen+1)] = EliteProcedure(OffspringPop,Population,BestSol(gen),GenParam);
            %save('temp','gen','Population','BestSol')
            fprintf('\n Generation %d Best perfprmance %4.3f with DV %f',gen, BestSol(gen+1).Mremoved,sum(BestSol(gen+1).DV))
        end
        save("Iridium33MARun"+num2str(RunIDX)+LowLevelType,"BestSol","wDV","wM","LowLevelType","GenParam")
        fprintf('\nRun %d LowLevelType %s end ', RunIDX, LowLevelType)
    end
    Lstart = 2;
end

