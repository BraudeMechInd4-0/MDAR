function Population = CalcRankAndDistance(Population)
%CALCRANKANDDISTANCE  Assign NSGA-II non-domination rank and crowding distance.
%   Population = CALCRANKANDDISTANCE(Population) forms the two objectives (negated
%   removed mass and total DeltaV), min-max normalizes them, performs fast
%   non-dominated sorting (NDSort), and computes the per-front crowding distance
%   (CrowdingDistance).
%
%   Inputs:
%     Population - struct array with fields .Mremoved and .DV
%   Outputs:
%     Population - same array with added fields .Rank (front number) and
%                  .CrowdDis (crowding distance)
F1 = -1*[Population(:).Mremoved]';
F2 = sum([Population(:).DV],1)';
FitnessValues = [F1,F2];

maxF = max(FitnessValues,[],1);
minF = min(FitnessValues,[],1);
if all(maxF~=minF)
    FitnessValues = (FitnessValues - minF)./(maxF - minF);
end

[FrontNo,~] = NDSort(FitnessValues,length(Population));
Distance = CrowdingDistance(FitnessValues,FrontNo);

for i = 1 : length(Population)
    Population(i).Rank = FrontNo(i);
    Population(i).CrowdDis = Distance(i);
end
end