function [offspring1, offspring2] = SBX(parent1,parent2,CrossOverRate,Limits)
%SBX  Simulated Binary Crossover for real-valued vectors.
%   [offspring1,offspring2] = SBX(parent1, parent2, CrossOverRate, Limits) applies
%   simulated binary crossover to a pair of real-valued vectors (here the maneuver
%   wait times), producing children clamped to the given bounds.
%
%   Inputs:
%     parent1, parent2 - real-valued vectors of equal length
%     CrossOverRate    - probability of applying SBX (0..1)
%     Limits           - [yl yu] lower/upper bounds applied to every gene
%   Outputs:
%     offspring1, offspring2 - SBX children, clamped to [yl,yu]
%
%   Note: the distribution index mu (eta_c, default 20) controls the spread;
%   larger values keep children closer to the parents. If crossover is not applied
%   the parents are returned unchanged.

mu = 20;
offspring1 = parent1;
offspring2 = parent2;
if rand(1) <= CrossOverRate

    yl = Limits(1);
    yu = Limits(2);
    for j = 1 : length(parent2)
        par1 = parent1(j);
        par2 = parent2(j);
        rnd = rand(1);
        if rnd <= 0.5
            if abs(par1 - par2) > 0.000001
                if par2 > par1
                    y2 = par2;
                    y1 = par1;
                else
                    y2 = par1;
                    y1 = par2;
                end
                if (y1 - yl) > (yu - y2)
                    beta = 1 + (2*(yu - y2)/(y2 - y1));
                else
                    beta = 1 + (2*(y1 - yl)/(y2 - y1));
                end
                expp = mu + 1;
                beta = 1/beta;
                alpha = 2 - beta^expp;
                rnd = rand(1);
                if rnd <= 1/alpha
                    alpha = alpha*rnd;
                    expp = 1/(mu + 1);
                    betaq = alpha^expp;
                else
                    alpha = alpha*rnd;
                    alpha = 1/(2 - alpha);
                    expp = 1/(mu + 1);
                    betaq = alpha^expp;

                end
                child1(j) = 0.5*((y1 + y2) - betaq*(y2 - y1));
                child2(j) = 0.5*((y1 + y2) + betaq*(y2 - y1));
            else
                betaq = 1;
                y1 = par1;
                y2 = par2;
                child1(j) = 0.5*((y1 + y2) - betaq*(y2 - y1));
                child2(j) = 0.5*((y1 + y2) + betaq*(y2 - y1));
            end
            if child1(j) < yl
                child1(j) = yl;
            end
            if child2(j) < yl
                child2(j) = yl;
            end
            if child1(j) > yu
                child1(j) = yu;
            end
            if child2(j) > yu
                child2(j) = yu;
            end
        else
            child1(j) = par1;
            child2(j) = par2;
        end
        if ~isreal(child1) || ~isreal(child2)
            disp('error in ga rep sbx');
            if ~isreal(child1)
                child1(j) = rand(1)*(yu - yl) + yl;
            end
            if ~isreal(child2)
                child2(j) = rand(1)*(yu - yl) + yl;
            end
        end
    end
    offspring1 = child1;
    offspring2 = child2;

else
end