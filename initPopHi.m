function Population = initPopHi(PopSizeHi,NumofDeb,mode,NtoRemove,T1Max)
%INITPOPHI  Initialize the upper-level (high-level) population.
%   Population = INITPOPHI(PopSizeHi, NumofDeb, mode, NtoRemove, T1Max) creates a
%   random initial population of removal plans (a debris visiting order plus the
%   maneuver wait times).
%
%   Inputs:
%     PopSizeHi - number of individuals in the population
%     NumofDeb  - total number of debris objects available
%     mode      - wait-time initialization: 1 -> random wait times in [0,T1Max];
%                 any other value leaves .WaitUntilManuver empty (set by the caller)
%     NtoRemove - number of debris selected per individual
%     T1Max     - maximum wait time between maneuvers [s]
%   Outputs:
%     Population - struct array with fields:
%                    .Order            - permutation of NtoRemove debris IDs
%                    .WaitUntilManuver - (NtoRemove-1)x1 wait times (if mode==1)
Population(PopSizeHi).Order = [];
Population(PopSizeHi).WaitUntilManuver = [];
for i = 1 : PopSizeHi
    Population(i).Order = randperm(NumofDeb,NtoRemove);
    if mode == 1
        Population(i).WaitUntilManuver = rand(NtoRemove-1,1)*T1Max;
    end
end