function [a,e,i,O,w,f] = kep_elements(r,v,mu)
%KEP_ELEMENTS  Classical Keplerian orbital elements from position and velocity.
%   [a,e,i,O,w,f] = KEP_ELEMENTS(r, v, mu) converts an ECI state to classical
%   orbital elements. Operates on single 1x3 position/velocity vectors.
%
%   Inputs:
%     r  - position vector [km] (1x3)
%     v  - velocity vector [km/s] (1x3)
%     mu - gravitational parameter [km^3/s^2]
%   Outputs:
%     a - semi-major axis [km]
%     e - eccentricity [-]
%     i - inclination [rad]
%     O - RAAN (right ascension of ascending node) [rad]
%     w - argument of periapsis [rad]
%     f - true anomaly [rad]

%% known vectors:
%ECI reference frame
iv = [1 0 0];
%jv = [0 1 0];
kv = [0 0 1];

norm_v = norm(v);
norm_r = norm(r);
%E = norm_v^2/2-mu/norm_r;
%h
h = cross(r,v);
norm_h = norm(h);
%n
n = cross(kv,h);
norm_n = norm(n);

%% a
a = (2/norm_r -norm_v^2/mu)^-1;

%% e
%vector
ev = cross(v,h)/mu-r/norm_r;
norm_ev = norm(ev);
e = norm(ev);

%% i
i = acos(dot(kv,h)/norm_h);

%% O
O = acos(dot(iv,n)/norm_n);
if n(2) < 0
    O = 2*pi - O;
end

%% w
w = acos(dot(n,ev)/(norm_n*norm_ev));
if ev(3) < 0
    w = 2*pi - w;
end

%% f
f = acos(dot(ev,r)/(norm_ev*norm_r));
if dot(r,v) < 0
    f = 2*pi - f;
end


if nargout == 1
    a = [a,e,i,O,w,f];
end


