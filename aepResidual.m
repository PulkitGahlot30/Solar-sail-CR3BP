function F = aepResidual(vars, mu, beta, alpha, delta)
% AEPRESIDUAL  Residual of the equilibrium conditions for the solar sail
% CR3BP, restricted to planar motion (z = 0), used with fsolve to find
% Artificial Equilibrium Points (AEPs) near the triangular points.
%
% At equilibrium, velocity and acceleration vanish in the synodic frame:
%   Omega_x + a_sail_x = 0
%   Omega_y + a_sail_y = 0
%
% INPUT vars = [x; y]  (planar AEP candidate)
% OUTPUT F   = [F1; F2] residuals (should be ~0 at the true AEP)

    x = vars(1);
    y = vars(2);
    z = 0;

    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

    Omega_x = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    Omega_y = y - (1-mu)*y/r1^3 - mu*y/r2^3;

    a_sail = solarSailAccel(x, y, z, mu, beta, alpha, delta);

    F = [Omega_x + a_sail(1);
         Omega_y + a_sail(2)];
end
