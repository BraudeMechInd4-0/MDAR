
%==========================================================================
% Bi-level Optimization Framework for MDAR (Proof-of-Concept)
%
% This script runs a high-level genetic optimization loop that evaluates
% debris removal missions (MDAR) using different low-level optimizers:
%   - "Rand" : Random search
%   - "GA"   : Genetic Algorithm
%   - "CMA"  : Covariance Matrix Adaptation
%
% Workflow:
%   1) Load models, constants, and Iridium-33 debris data.
%   2) Initialize high-level population and evaluate (calls low-level solver).
%   3) Iterate GA operations: selection, crossover, mutation, evaluation.
%   4) Track/keep elite solution each generation.
%   5) Save best solutions per run and low-level type.
%
% Notes:
%   - Requires subfolders: ParseGP, Satellites, vallado, models.
%   - Uses pre-saved MATLAB data files: models/Iridium33Model.mat, Iridium33.mat.
%   - Weights wM and wDV control fitness tradeoff (removed mass vs ΔV).
%==========================================================================

clc; clear; close all                      % Reset environment
rng('shuffle')                             % Randomize RNG seed

%% Setup and data loading
fprintf('loading...');
addpath(genpath("ParseGP"))                % Add parsing utilities
addpath(genpath("Satellites"))             % Add satellite-related utilities
addpath(genpath("vallado"))                % Add orbital mechanics utilities
addpath(genpath("models"))                 % Add models folder

load("models\Iridium33Model.mat")          % Load Iridium-33 model parameters
Cd = 2.2;                                  % Drag coefficient
mu1 = 398600.5;                            % Earth GM [km^3/s^2]
Re = 6378.137;                             % Earth mean radius [km]
J = [0, 0.01082635854, ...                 % Zonal harmonics (J2..J6-like set)
    -0.2532435346e-5, -0.1619331205e-5, ...
    -0.2277161016e-6 , 0.5396484906e-6];

load('Iridium33.mat');                     % Load debris set and epochs
fprintf("done!\n")

%% High-level GA parameters
PopSizeHi = 30;                            % Population size (high level)
MaxGenHi  = 50;                            % Number of generations (high level)
mc        = 500;                           % Monte Carlo samples (low level)
Ac        = 9e-6;                          % Cross-sectional area [km^2]
NtoRemove = 10;                            % Number of debris to remove in a plan
T1Max     = 2*60*60*24*30;                 % Max Waiting time between manuvers [s] (≈60 days)

GenParam.pm        = 0.4;                  % Mutation probability
GenParam.pc        = 0.8;                  % Crossover probability
wM                 = 1;                    % Weight for removed mass metric
wDV                = 8;                    % Weight for total ΔV metric
GenParam.PopSize   = PopSizeHi;            
GenParam.T1max     = T1Max;
GenParam.NtoRemove = NtoRemove;
GenParam.NDeb      = length(m);            % Number of debris items
GenParam.EliteMode = "H";                  % High-level elite handling mode

[~, ILastUpdated]  = min(days(datetime("now") - epoch)); % Find most recent epoch index
RelativeEpoch      = seconds(epoch(ILastUpdated) - epoch);% Time offset vector [s]

Runs    = 3;                               % Number of outer runs
Lstart  = 1;                               % Start Ltype (low-level solver types) index

%% Outer runs over low-level solver types
for RunIDX = 1:Runs                         
    for Ltype = Lstart : 3                  % Iterate low-level type (1..3)
        % Map Ltype to LowLevelType string
        if Ltype == 1
            LowLevelType = "Rand";          % Random search low-level
        elseif Ltype == 2
            LowLevelType = "GA";            % GA low-level
        else
            LowLevelType = "CMA";           % CMA-ES low-level
        end

        fprintf('\nRun %d LowLeverType %s Start ... ', RunIDX, LowLevelType)

        %---------- Initialization and first evaluation ----------
        Population = initPopHi(PopSizeHi, length(m), NtoRemove, T1Max);   % Create initial population
        Population = EvaluateHi(Population, model, m, RelativeEpoch, ...  % Evaluate with low-level optimizer
            mc, Ac, J, Cd, Re, mu1, wM, wDV, LowLevelType);

        gen = 0;                                                          % Generation counter
        [~, idxBest] = max([Population(:).fit]);                          % Get best individual index
        BestSol = Population(idxBest);                                    % Store best solution
        Ifmax = 1;                                                        % Flag: maximize fitness

        fprintf('\n Generation %d Best perfprmance %4.3f with DV %f', ...
            gen, BestSol(1).Mremoved, sum(BestSol(1).DV))             % Report baseline

        %---------- Evolution loop ----------
        for gen = 1 : MaxGenHi
            PopulationSel = tournamentSelection(Population, Ifmax);       % Select parents
            PopulationCO  = CrossOver(PopulationSel, GenParam);           % Apply crossover
            OffspringPop  = Mutation(PopulationCO, GenParam);             % Apply mutation

            % Evaluate offspring using the chosen low-level optimizer
            OffspringPop  = EvaluateHi(OffspringPop, model, m, RelativeEpoch, ...
                mc, Ac, J, Cd, Re, mu1, wM, wDV, LowLevelType);

            % Elite procedure: merge offspring + previous population; keep best
            [Population, BestSol(gen+1)] = EliteProcedure(OffspringPop, ...
                Population, ...
                BestSol(gen), ...
                GenParam);

            % Progress report
            fprintf('\n Generation %d Best perfprmance %4.3f with DV %f', ...
                gen, BestSol(gen+1).Mremoved, sum(BestSol(gen+1).DV))
        end

        %---------- Save results for this run/type ----------
        save("Run" + num2str(RunIDX) + LowLevelType, ...
            "BestSol", "wDV", "wM", "LowLevelType", "GenParam")

        fprintf('\nRun %d LowLevelType %s end ', RunIDX, LowLevelType)
    end

    Lstart = 1;                                  % Reset Ltype start for next run
end
