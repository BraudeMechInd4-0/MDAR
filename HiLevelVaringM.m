%HILEVELVARINGM  Weighted-sum (scalarization) baseline for the multi-objective study.
%   Sweeps a set of objective weight vectors and, for each, runs the bi-level GA with
%   the mass-varying evaluator (EvaluateHiVarM) to trace a front by scalarization — a
%   comparison baseline for the NSGA-II results (paper 2).
%
%   Uses the J6+Drag+SRP+Moon force model and the Thompson lower level (prtlambertT).
%   For each weight vector: initialize the population, evaluate, then iterate
%   selection/crossover/mutation/evaluation/elitism; save per-weight, per-run results.
%
%   Requires: ParseGP, Satellites, fundamentals-of-astrodynamics, models. Loads the
%   precomputed <cloud><forceModel>Model.mat and the <cloud>.mat catalog.
%   Key parameters (below): PopSizeHi, MaxGenHi, NtoRemove, T1Max, weight set WM.
clc;clear;close all
rng('shuffle')
%% Weighted-sum scalarization sweep
% General Parameters + propagating the Satellite locaions
parfevalOnAll(@warning, 0, 'off', 'all');
fprintf('loading...');
addpath(genpath("ParseGP"))
addpath(genpath("Satellites"))
addpath(genpath("fundamentals-of-astrodynamics"))
addpath(genpath("models"))
forceModel = "J6DSM";
for C = 1:1:2
     
    if C == 1
        DebrisCloud = "Iridium33"; % Iridium33 OR hazardous50 OR SyntheticCloud
	WM = [-0.01,-0.0231,0.0123];
    else
        DebrisCloud = "Iridium33"; % Iridium33 OR hazardous50 OR SyntheticCloud
	WM = [-0.08,-0.2,0.12];
    end
    %load("models\Iridium33Model.mat")
    %load("models\hazardous50Model.mat")
    load("models\"+ DebrisCloud +forceModel+"Model.mat")


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
    %load('Iridium33.mat');
    load(DebrisCloud + ".mat")
    %save Satellites\IriduimList.mat m A epoch r v
    fprintf("done!\n")

    %RunIDX = 2;
    PopSizeHi = 30;
    MaxGenHi = 50;
    mc = 500;%m(10);
    Ac = 9e-6;%A(10);


    if strcmp(DebrisCloud,"SyntheticCloud")
        RelativeEpoch = zeros(size(m));

    else

        [~,ILastUpdated] = min(days(datetime("now") - epoch));
        RelativeEpoch = seconds(epoch(ILastUpdated) - epoch);
    end
    NtoRemove = 10;
    T1Max =  2*60*60*24*30;
    GenParam.pm = 0.4;
    GenParam.pc = 0.8;
    %wM = 1;
    wDV = 1;
    elitePer = 1;
    GenParam.PopSize = PopSizeHi;
    GenParam.T1max = T1Max;
    GenParam.NtoRemove = NtoRemove;
    GenParam.NDeb = length(m);
    GenParam.EliteMode = "H";
    Data(MaxGenHi+1).Pop = [];
    Data(MaxGenHi+1).ND = [];

    Data(MaxGenHi+1).NDValue = [];
    Runs = 7;
    Isolver = 1;
    LowLevelType = "Rand";
    mode = 1;Lstart = 2;
	Ifmax = 0;
    for RunIDX = 1:Runs
	
	%BestSol = struct;
	SolverString = "T";
    Solver = 1;
        for idxW = 1:3
        Data = struct;
	    wM = WM(idxW);
            fprintf('\nRun %d LowLeverType %s Start ... ', RunIDX, LowLevelType)
            tic
            Population = initPopHi(PopSizeHi,length(m),mode,NtoRemove,T1Max);
            %
            Population = EvaluateHiVarM(Population,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType,Solver);
            %toc
            gen = 0;
           
            Fit = [Population(:).fit];
            [M, I] = min(Fit);
	        BestSol(idxW) = Population(I);
            
            fprintf('\n Generation %d Best DV  %4.3f Best M %4.3f Best Fit %5f',gen, sum(BestSol(idxW).DV),BestSol(idxW).Mremoved,BestSol(idxW).fit)
            Data(gen+1).Pop = Population;
            Data(gen+1).BestSol = BestSol(idxW);
            
            for gen = 1 : MaxGenHi
                
                PopulationSel = tournamentSelection(Population,Ifmax);
                PopulationCO = CrossOver(PopulationSel,GenParam);
                OffspringPop = Mutation(PopulationCO,GenParam);
                OffspringPop = EvaluateHiVarM(OffspringPop,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType,Solver);

                [Population, BestSol(idxW)] = EliteProcedure(OffspringPop,Population,BestSol(idxW),GenParam);
                
            
                
                fprintf('\n Generation %d Best DV  %4.3f Best M %4.3f Best Fit %5f',gen, sum(BestSol(idxW).DV),BestSol(idxW).Mremoved,BestSol(idxW).fit)
                Data(gen+1).Pop = Population;
                Data(gen+1).BestSol = BestSol(idxW);
                
                
            end
            runTime = toc;

            %save("Iridium33VMRun"+num2str(RunIDX)+LowLevelType+"Nremoved"+NtoRemove,"BestSol","wDV","wM","LowLevelType","GenParam")
            save(DebrisCloud + forceModel + "VM_SO_W_"+num2str(idxW)+"_Run"+num2str(RunIDX)+LowLevelType+"Nremoved"+NtoRemove,"Data","GenParam","runTime")
            fprintf('\nRun %d LowLevelType %s end , Time = %e', RunIDX, LowLevelType,runTime)
        end
        clear BestSol
    end

end