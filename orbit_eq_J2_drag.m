function [dr] = orbit_eq_J2_drag(t,r,mu1,CD,A,m,Re,J2)
%ORBIT_EQ_J2_DRAG  Equations of motion: J2 zonal harmonic + atmospheric drag.
%   dr = ORBIT_EQ_J2_DRAG(t, r, mu1, CD, A, m, Re, J2) returns the state derivative
%   for a spacecraft under Earth gravity with the J2 zonal harmonic and drag.
%
%   Inputs:
%     t   - time [s] (ODE-solver argument)
%     r   - state [rx ry rz vx vy vz] (6xN column-wise, or Nx6 row-wise)
%     mu1 - Earth gravitational parameter [km^3/s^2]
%     CD  - drag coefficient
%     A   - cross-sectional area [km^2]
%     m   - spacecraft mass [kg]
%     Re  - Earth radius [km]
%     J2  - J2 zonal harmonic coefficient
%   Outputs:
%     dr  - state derivative, same shape as r


% d^2r/dt = -mu*r/norm(r);

[n,M] = size(r);
dr = zeros (size(r));


if n == 6

    normr = sqrt(sum(r(1:3,:).^2,1));

    ad = drag_accel(r',CD,A,m,Re)';

    aj2 = -3/2*J2*(mu1./(normr.^2)).*(Re./normr).^2.*[(1-5*(r(3,:)./normr).^2).*r(1,:)./normr, (1-5*(r(3,:)./normr).^2).*r(2,:)./normr, (3-5*(r(3,:)./normr).^2).*r(3,:)./normr]';
    % ad = zeros(size(aj2));
    le = normr.^3;

    dr(1,:) = r(4,:);
    dr(2,:) = r(5,:);
    dr(3,:) = r(6,:);
    dr(4,:) = -mu1*r(1,:)./le + ad(1,:)+aj2(1,:);
    dr(5,:) = -mu1*r(2,:)./le + ad(2,:)+aj2(2,:);
    dr(6,:) = -mu1*r(3,:)./le + ad(3,:)+aj2(3,:);
elseif M == 6

    normr = sqrt(sum(r(:,1:3).^2,2));

    ad = drag_accel(r,CD,A,m,Re);
    aj2 = -3/2*J2*(mu1./(normr.^2)).*(Re./normr).^2.*[(1-5*(r(:,3)./normr).^2).*r(:,1)./normr, (1-5*(r(:,3)./normr).^2).*r(:,2)./normr, (3-5*(r(:,3)./normr).^2).*r(:,3)./normr];
    % ad = zeros(size(aj2));
    le =normr.^3;

    dr(:,1) = r(:,4);
    dr(:,2) = r(:,5);
    dr(:,3) = r(:,6);
    dr(:,4) = -mu1*r(:,1)./le + ad(:,1)+aj2(:,1);
    dr(:,5) = -mu1*r(:,2)./le + ad(:,2)+aj2(:,2);
    dr(:,6) = -mu1*r(:,3)./le + ad(:,3)+aj2(:,3);
end






