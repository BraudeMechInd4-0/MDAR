function offspring = PolyMutation(parent,MutationRate,Limits)
%POLYMUTATION Apply polynomial mutation to real-valued genes.
% Uses polynomial probability distribution to perturb genes near parent values.
%
% Inputs:
%   parent       - real-valued vector (genes to mutate)
%   MutationRate - probability of applying mutation [0,1]
%   Limits       - [lower, upper] bounds for genes
% Outputs:
%   offspring    - mutated vector (same size as parent)
%
% Notes:
%   - Distribution index mu=20 controls perturbation spread
%   - Based on Deb & Agrawal (1995) polynomial mutation operator

Nl = length(parent);
mu = 20;
offspring = parent;
if rand(1) <= MutationRate
    for j = 1 : Nl
        yl = Limits(1);
        yu = Limits(2);
        y = parent(j);
        if y > yl
            if (y - yl) < (yu - y)
                delta = (y - yl)/(yu - yl);
            else
                delta = (yu - y)/(yu - yl);
            end
            rnd = rand(1);
            indi = 1/(mu + 1);

            if rnd <= 0.5
                u = 1 - delta;
                val = 2*rnd + (1-2*rnd)*(u^(mu+1));
                deltaq = val^indi - 1;
            else
                u = 1 - delta;
                val = 2*(1 - rnd) + 2*(rnd - 0.5)*(u^(mu + 1));
                deltaq = 1 - val^indi;
            end
            y = y + deltaq * (yu - yl);
            if y > yu
                y = yu;
            end
            if y < yl
                y = yl;
            end
            
             offspring(j) = y;
            
        else
            
              offspring(j) = rand(1)*(yu - yl) + yl;
            
        end
        if ~isreal(offspring(j))
            warning('error in ga rep polmutation');
            offspring(j) = rand(1)*(yu - yl) + yl;
        end
    end
end