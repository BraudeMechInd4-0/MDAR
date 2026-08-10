function PopulationAfterCross = CrossOver(Population,GenParam)
%CROSSOVER  Recombine a population: SBX on wait times, OX on visiting order.
%   PopulationAfterCross = CROSSOVER(Population, GenParam) recombines individuals
%   pairwise ((1,2),(3,4),...): simulated binary crossover (SBX) on the maneuver
%   wait times and ordered crossover (OX) on the debris visiting sequence,
%   controlled by GenParam.pc.
%
%   Inputs:
%     Population - struct array with fields .WaitUntilManuver, .Order
%     GenParam   - struct with fields:
%                    .pc      - crossover probability (0..1)
%                    .NDeb    - number of debris labels (permutation repair)
%                    .PopSize - population size
%                    .T1max   - upper bound for wait times (lower bound = 1)
%   Outputs:
%     PopulationAfterCross - population after crossover
pc = GenParam.pc;
NDeb = GenParam.NDeb;
PopulationAfterCross = Population;
PopSize = GenParam.PopSize;
T1limits = [1 GenParam.T1max];
for i = 1 : 2 :PopSize
    if rand() > pc
        [PopulationAfterCross(i).WaitUntilManuver,PopulationAfterCross(i+1).WaitUntilManuver] = SBX(Population(i).WaitUntilManuver,Population(i+1).WaitUntilManuver,1,T1limits);
        [PopulationAfterCross(i).Order,PopulationAfterCross(i+1).Order] = Crossover_Ordered_Operator(Population(i).Order,Population(i+1).Order,NDeb);

    else
        PopulationAfterCross(i).WaitUntilManuver = Population(i).WaitUntilManuver;
        PopulationAfterCross(i+1).WaitUntilManuver = Population(i+1).WaitUntilManuver;
        PopulationAfterCross(i).Order = Population(i).Order;
        PopulationAfterCross(i+1).Order = Population(i+1).Order;
    end

end