function OrderAfter = ExchangeMutation(Order,pm,NDeb)
%EXCHANGEMUTATION Apply exchange or replacement mutation to permutation.
% Randomly swaps two positions (30% prob) or replaces one gene with unused label.
%
% Inputs:
%   Order   - permutation vector (debris visit sequence)
%   pm      - mutation probability [0,1]
%   NDeb    - maximum debris label (for sampling unused labels)
% Outputs:
%   OrderAfter - mutated permutation (same length as Order)
%
% Strategy:
%   With prob pm:
%     30% chance: Exchange mutation (swap two random positions)
%     70% chance: Replacement mutation (replace one gene with unused label)

if rand() < pm
    if rand() < 0.3
       %do exchange
       S = randi(length(Order),1,2);
       OrderAfter = Order;
       OrderAfter(S(1)) = Order(S(2));
       OrderAfter(S(2)) = Order(S(1));
    else
        %do replace
        S = randi(length(Order),1,1);
        OrderAfter = Order;
        s = randi(NDeb,1,1);
        while any(s == Order)
            s = randi(NDeb,1,1);
        end
        OrderAfter(S(1)) = s;
        
    end
else
OrderAfter = Order;
end
end