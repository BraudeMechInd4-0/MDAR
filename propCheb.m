function xout = propCheb(t,model)
%This function propagates a satellite using a pre-computed chebyshev
%polynomail.
%Usage:
% xout = propCheb(t,model)
%
%Where:
%   xout - is the satellite state, a size 6 vector of position and velocity
%   t    - is the time we wish to propagate to
%   model - a structure containing the model for the satellite
%   model.t - the times in which each polynomial starts
%   model.bi - the polynomial coeeficients acc. to woolands and bai's works

%First find the right polynomials 
ind = find (t > model.ts,1, 'last');
if isempty(ind) || ind == length(model.ts)
    error("point not in model");
end
tstart = model.ts(ind);
tend   = model.ts(ind+1);
om2 = (tend-tstart)/2;
om1 = (tend+tstart)/2;
bi(:,:) = model.bi(ind,:,:);
N = size(bi,1)-1;
tau = (t-om1)/om2;
T = cos((0:N)'*acos(tau))';

xout = T*bi;
