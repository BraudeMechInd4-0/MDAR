
function PopulationAfterMut = Mutation(Population, GenParam)
%MUTATION Apply polynomial mutation to wait times and exchange mutation to orders.
% Inputs:
%   Population - struct array with fields .WaitUntilManuver, .Order
%   GenParam   - struct with fields:
%                .pm     -> mutation probability (0..1)
%                .PopSize-> population size
%                .T1max  -> upper bound for wait times
%                .NDeb   -> number of debris labels (for repair range)
% Output:
%   PopulationAfterMut - population after mutation

    pm   = GenParam.pm;                 % mutation probability
    PopSize = GenParam.PopSize;         % population size
    T1limits = [1 GenParam.T1max];      % bounds for wait times (polynomial mutation)
    NDeb = GenParam.NDeb;               % label max (for exchange mutation bounds)

    PopulationAfterMut = Population;    % initialize output

    % Mutate each individual independently
    for i = 1:PopSize
        if rand() < pm
            % Apply polynomial mutation on continuous wait times
            PopulationAfterMut(i).WaitUntilManuver = ...
                PolyMutation(Population(i).WaitUntilManuver, 1, T1limits);

            % Apply exchange mutation on permutation of debris order
            PopulationAfterMut(i).Order = ...
                ExchangeMutation(Population(i).Order, 1, NDeb);
        else
            % No mutation: carry over parent as-is
            PopulationAfterMut(i).WaitUntilManuver = Population(i).WaitUntilManuver;
            PopulationAfterMut(i).Order            = Population(i).Order;
        end
    end
end
