function PopulationAfterMut = Mutation(Population,GenParam)
%MUTATION  Mutate a population: polynomial mutation on waits, exchange on order.
%   PopulationAfterMut = MUTATION(Population, GenParam) applies polynomial mutation
%   to the maneuver wait times and exchange (swap) mutation to the debris visiting
%   sequence, controlled by GenParam.pm.
%
%   Inputs:
%     Population - struct array with fields .WaitUntilManuver, .Order
%     GenParam   - struct with fields:
%                    .pm      - mutation probability (0..1)
%                    .PopSize - population size
%                    .T1max   - upper bound for wait times
%                    .NDeb    - number of debris labels
%   Outputs:
%     PopulationAfterMut - population after mutation
pm = GenParam.pm;
PopulationAfterMut = Population;
PopSize = GenParam.PopSize;
T1limits = [1 GenParam.T1max];
NDeb = GenParam.NDeb;
for i = 1 :PopSize
    if rand() > pm
        PopulationAfterMut(i).WaitUntilManuver = PolyMutation(Population(i).WaitUntilManuver,1,T1limits);
        PopulationAfterMut(i).Order = ExchangeMutation(Population(i).Order,1,NDeb);

    else
        PopulationAfterMut(i).WaitUntilManuver = Population(i).WaitUntilManuver;
        PopulationAfterMut(i).Order = Population(i).Order;
    end

end

