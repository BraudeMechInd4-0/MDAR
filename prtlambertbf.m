function [v1t, v2t, errorb] = prtlambertbf(forcemodel,odesolver, r1, r2, v1, dm, de, nrev, dtsec,mu1, Tolr, solveroptions )
%
% solve the pertubed lamberts problem. 
%Usage:
%    [v1t, v2t, errorb] = prtlambertbf(forcemodel,odesolver, r1, r2, v1, dm, de, nrev, dtsec,mu1 )
%
% Inputs:
%   forcemodel - a function handle to the force model, should be defined
%   with the parameters, examle:
%       @(t,x)orbit_eq_J2_drag(t,x,mu,Cd,At,mt,Re,J2)
%   odesolver - a function handle to the odesolver, should be something
%   like @ode45
%   r1 - initial position       km
%   r2 - end position           km
%   v1 - initial velocity       km/s
%   dm - direction of motion ('L', 'S')
%   de - orbital energy ('L', 'H')
%   nrev - number of revs to complete 
%   dtsec - time between r1 and r2          s
%
% Outputs:
%    v1t - ijk transfer velocity vector               km/s
%    v2t - ijk transfer velocity vector               km/s
maxIter = 20;
deltar = inf;
prop_r2 = r2;
iterCounter = 1;
flag = 0;
while norm(deltar) >= Tolr && iterCounter < maxIter
    [v1t,v2t,errorb] = lambertb ( r1,prop_r2,v1, dm,de,nrev, dtsec,mu1);
    if ~isreal(v1t) || ~isreal(v2t) 
        %warning("could not solve, bad parameters"); 
        break
    end
    if  hitearth ( 100, r1, v1t, r2, v2t, nrev )
        flag = 1;
        errorb = 'Hit earth';
        break
    end
    [~,x2]=odesolver(forcemodel,[0 dtsec],[r1 v1t],solveroptions);
    r_new = x2(end,1:3);
    v2t = x2(end,4:6);
    deltar = r_new - r2;
    prop_r2 = prop_r2-deltar;
    iterCounter = iterCounter + 1;
end

if iterCounter == maxIter || flag
    v1t = v1t + 1i;
    v2t = v2t + 1i;
end