function Population = initPopHi(PopSizeHi, NumofDeb, NtoRemove, T1Max)
%INITPOPHI Initialize high-level GA population
% Inputs:
%   PopSizeHi  - number of individuals in the population
%   NumofDeb   - total number of debris items available
%   NtoRemove  - number of debris to select per individual
%   T1Max      - maximum waiting time between maneuvers [s]
% Output:
%   Population - struct array with fields:
%                .Order            -> permutation of selected debris IDs
%                .WaitUntilManuver -> random wait times between maneuvers

    %--- Preallocate struct array for efficiency ---
    Population(PopSizeHi).Order = [];
    Population(PopSizeHi).WaitUntilManuver = [];

    %--- Generate each individual ---
    for i = 1:PopSizeHi
        Population(i).Order = randperm(NumofDeb, NtoRemove);      % Random debris order
        Population(i).WaitUntilManuver = rand(NtoRemove-1, 1) * T1Max; % Random wait times
    end
end
