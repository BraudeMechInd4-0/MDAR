%MDAR_NSGAII_MAIN  Multi-objective bi-level MDAR optimization (paper 2, NSGA-II).
%   NSGA-II upper level over the debris visiting sequence and maneuver timing, with
%   two objectives: maximize removed mass and minimize total DeltaV. The lower level
%   (Thompson perturbed-Lambert) resolves each transfer under the J6+Drag+SRP+Moon
%   force model. Runs the Iridium-33 and synthetic-cloud scenarios, repeated to pool
%   a reference Pareto front; saves per-run non-dominated sets and per-generation logs.
%
%   Requires: ParseGP, Satellites, fundamentals-of-astrodynamics, models.
%   Key parameters (below): PopSizeHi, MaxGenHi, NtoRemove, T1Max, Runs.
clc;clear;close all
rng('shuffle')                                              % randomize the RNG seed so successive runs are not identical
%% This program runs a NSGAII-based for Designing Multi-Debris Active Removal Missions
% General Parameters + propagating the Satellite locaions
% -------------------------------------------------------------------
% Implements the multi-objective bi-level framework (EMOBL) described
% in the paper: an NSGA-II upper level searches over the debris
% visiting sequence + maneuver timing (two objectives: maximize
% removed mass, minimize total DeltaV), while a lower-level solver
% (Thompson's perturbed-Lambert method, driven by random search)
% resolves the transfer parameters for each leg and reports DeltaV
% back to the upper level.
% -------------------------------------------------------------------
parfevalOnAll(@warning, 0, 'off', 'all');                   % silence warnings on all parallel-pool workers (e.g., Lambert-solver edge cases)
fprintf('loading...');
addpath(genpath("ParseGP"))                                 % upper-level GA/NSGA-II operators (init, crossover, mutation, sorting, etc.)
addpath(genpath("Satellites"))                               % satellite/debris propagation utilities
addpath(genpath("fundamentals-of-astrodynamics"))                                  % Vallado astrodynamics routines (Lambert solvers, ephemerides, etc.)
addpath(genpath("models"))                                   % precomputed/propagated debris-cloud + force-model data
forceModel = "J6DSM";                                        % force-model tag: J6 zonal harmonics + Drag + SRP + Moon (3rd body)
for C = 2:-1:1                                                % loop over the two test scenarios: C=2 -> Iridium33 (real cloud), C=1 -> SyntheticCloud
    if C == 1
        DebrisCloud = "SyntheticCloud"; % Iridium33 OR hazardous50 OR SyntheticCloud
    else
        DebrisCloud = "Iridium33"; % Iridium33 OR hazardous50 OR SyntheticCloud
    end
    load("models\"+ DebrisCloud +forceModel+"Model.mat")      % load precomputed propagation model for this cloud + force model (provides "model", and, based on later usage, "Solver")
    Cd = 2.2;                                                 % drag coefficient (cannonball drag model, Appendix A)
    mu1 = 398600.5;                                           % Earth's gravitational parameter, km^3/s^2
    Re = 6378.137;                                            % Earth's radius, km
    J = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, -0.2277161016e-6 ,0.5396484906e-6]; % zonal harmonic coefficients [J1..J6] used to build the a_z acceleration up to J6
    %% generate list of objects and parameters: now it is for the chaser and one target
    load(DebrisCloud + ".mat")                                % load debris-cloud catalog (orbital elements, masses "m", areas "A", TLE "epoch", etc.)
    fprintf("done!\n")
    PopSizeHi = 30;                                           % upper-level (NSGA-II) population size
    MaxGenHi = 50;                                            % upper-level number of generations
    mc = 500;%m(10);                                          % chaser (servicing spacecraft) mass, kg
    Ac = 9e-6;%A(10);                                          % chaser cross-sectional area, km^2 (= 9 m^2)
    if strcmp(DebrisCloud,"SyntheticCloud")
        RelativeEpoch = zeros(size(m));                       % synthetic cloud: all objects share the same reference epoch
    else
        [~,ILastUpdated] = min(days(datetime("now") - epoch));     % find the catalog object whose TLE epoch is closest to "now"
        RelativeEpoch = seconds(epoch(ILastUpdated) - epoch);      % express every object's epoch relative to that most-recently-updated epoch
    end
    NtoRemove = 10;                                           % number of debris objects to select and remove (m in the paper)
    T1Max =  2*60*60*24*30;                                   % maximum allowed inter-maneuver wait time, seconds (2 months, per the upper-level encoding bound)
    GenParam.pm = 0.4;                                         % mutation probability
    GenParam.pc = 0.8;                                         % crossover probability
    wM = 1; %not used here                                     % legacy single-objective mass weight from the scalarized formulation (Eq. 3); kept only for function-interface compatibility
    wDV = 8; %not used here                                    % legacy single-objective DeltaV weight; unused in this multi-objective (NSGA-II) version
    elitePer = 1;                                              % elitism parameter passed to the environmental-selection routine
    GenParam.PopSize = PopSizeHi;
    GenParam.T1max = T1Max;
    GenParam.NtoRemove = NtoRemove;
    GenParam.NDeb = length(m);                                 % number of candidate debris objects in the catalog (M in the paper)
    GenParam.EliteMode = "H";                                  % elitism/selection mode flag used by EliteFullSorting
    Data(MaxGenHi+1).Pop = [];                                 % preallocate per-generation logging struct (index g+1 stores generation g, so generations 0..MaxGenHi fit in 1..MaxGenHi+1)
    Data(MaxGenHi+1).ND = [];                                  % non-dominated individuals at each generation
    Data(MaxGenHi+1).NDValue = [];                             % their objective-vector values
    Runs = 7;                                                  % number of independent repeated runs per scenario (used later to pool/build the reference MO-Front)
    Isolver = 1;                                               % solver index (not referenced again in this script)
    LowLevelType = "Rand";                                     % lower-level trajectory search strategy: random search over Lambert-branch parameters
    mode = 1;Lstart = 2;                                       % upper-level encoding options passed to initPopHi (Lstart is not passed on to any call below, so it is unused in this script)
    for RunIDX = 1:Runs                                            % repeat the full NSGA-II optimization "Runs" times (for later statistical pooling into the reference front)
        SolverString = "T"; %T - Thompson                       % lower-level Lambert solver: Thompson's perturbed-Lambert method (ref. [22] in the paper)
        fprintf('\nRun %d LowLeverType %s Start ... ', RunIDX, LowLevelType)
        tic                                                     % start this run's timer
        Population = initPopHi(PopSizeHi,length(m),mode,NtoRemove,T1Max);   % randomly initialize the upper-level population (removal sequence S + timing Ts for each individual)
        Population = EvaluateHiVarM(Population,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType,Solver);   % evaluate each individual: run the lower-level solver on every leg to get DeltaV, mass removed, transfer/wait times, etc.
        gen = 0;                                                % generation counter (0 = initial population)
        F1 = -1*[Population(:).Mremoved]';                     % objective 1 (to minimize): negated removed mass, so that maximize-mass becomes minimize(-mass)
        F2 = sum([Population(:).DV],1)';                       % objective 2 (to minimize): total DeltaV summed over all maneuvers/legs
        F3 = sum([Population(:).TM],1)' + sum([Population(:).WaitUntilManuver],1)';    % total mission duration (transfer durations + wait times); tracked/logged alongside the two optimized objectives
        FitnessValues = [F1,F2];                                % combined fitness matrix used for non-dominated sorting
        [FrontNo,~] = NDSort(FitnessValues,1);                  % fast non-dominated sorting (PlatEMO-style); the trailing "1" restricts sorting to identifying only the first front
        NonDominated = Population(FrontNo==1);                  % current non-dominated ("elite") individuals
        NonDominatedValues = FitnessValues(FrontNo==1,:);       % their objective vectors
        fprintf('\n Generation %d Best DV  %4.3f Best M %4.3f Best Time %5f',gen, min(F2),abs(min((F1))),min(F3))   % progress printout for generation 0
        Data(gen+1).Pop = Population;                           % log full population for this generation
        Data(gen+1).ND = NonDominated;                          % log non-dominated set for this generation
        Data(gen+1).NDValue = NonDominatedValues;               % log non-dominated objective values for this generation
        for gen = 1 : MaxGenHi                                          % main NSGA-II generational loop
            Population = CalcRankAndDistance(Population);        % assign non-domination rank and crowding distance to each individual (NSGA-II selection criteria)
            PopulationSel = SelectionByRank(Population);          % build mating pool via rank-then-crowding tournament selection
            PopulationCO = CrossOver(PopulationSel,GenParam);      % recombine parents (ordered crossover on the sequence, SBX on the continuous timing) using GenParam.pc
            OffspringPop = Mutation(PopulationCO,GenParam);        % mutate offspring (exchange mutation on the sequence, polynomial mutation on the timing) using GenParam.pm
            OffspringPop = EvaluateHiVarM(OffspringPop,model,m,RelativeEpoch,mc,Ac,J,Cd,Re,mu1,wM,wDV,LowLevelType,Solver);   % evaluate offspring via the lower-level solver
            Population = EliteFullSorting(NonDominated,OffspringPop,PopSizeHi,elitePer);   % elitist environmental selection: combine previous elites + offspring, keep the best PopSizeHi by rank/crowding
            F1 = -1*[Population(:).Mremoved]';                    % recompute objective 1 (negated removed mass) for the new population
            F2 = sum([Population(:).DV],1)';                      % recompute objective 2 (total DeltaV)
            F3 = sum([Population(:).TM],1)' + sum([Population(:).WaitUntilManuver],1)';   % recompute total mission duration (logged only)
            N = Population(:).Nremoved;                            % NOTE: unlike F1/F2/F3 above, this is not wrapped in [...]; for a struct array, Population(:).Nremoved is a comma-separated list, so this assignment keeps only the first individual's Nremoved rather than the full vector - worth double-checking this is intended, since max(N) is used just below for reporting
            FitnessValues = [F1,F2];                               % combined fitness matrix for this generation
            [FrontNo,~] = NDSort(FitnessValues,1);                 % non-dominated sorting of the new population
            NonDominated = Population(FrontNo==1);                 % updated non-dominated set (elites carried into next generation)
            NonDominatedValues = FitnessValues(FrontNo==1,:);      % their objective vectors
            Data(gen+1).Pop = Population;                          % log this generation's population
            Data(gen+1).ND = NonDominated;                         % log this generation's non-dominated set
            Data(gen+1).NDValue = NonDominatedValues;              % log this generation's non-dominated objective values
            fprintf('\n Generation %d Best DV  %4.3f Best M %4.3f Best Time %5f with %d removed',gen, min(F2),abs(min((F1))),min(F3), max(N))   % progress printout for this generation
        end
        runTime = toc;                                            % stop this run's timer
        save(DebrisCloud + forceModel + "VM_MO_Three_Run"+num2str(RunIDX)+LowLevelType+"Nremoved"+NtoRemove,"NonDominated","Data","GenParam","runTime")   % save this run's results, filename encodes cloud, force model, run index, lower-level type, and NtoRemove
        fprintf('\nRun %d LowLevelType %s end , Time = %e', RunIDX, LowLevelType,runTime)
    end
end
