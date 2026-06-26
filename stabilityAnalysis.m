function [is_stable, eigenvalues, J] = stabilityAnalysis(aep_pos, mu, beta, alpha, delta)
% STABILITYANALYSIS  Linear stability of an Artificial Equilibrium Point
% (AEP) in the planar solar sail CR3BP, via the standard 4x4 state-space
% Jacobian of the linearized variational equations:
%
%   [dxi]     [ 0      0      1      0  ] [xi ]
%   [deta]  = [ 0      0      0      1  ] [eta]
%   [dxid]    [Uxx    Uxy     0      2  ] [xid]
%   [detad]   [Uxy    Uyy    -2      0  ] [etad]
%
% where Uxx, Uxy, Uyy are second partials of the TOTAL potential
% (gravity + sail) evaluated at the AEP. Because the sail acceleration
% is a nonlinear function of position (through r_s, n_hat, etc.), these
% second partials are computed via central finite differences rather
% than closed-form differentiation.
%
% INPUTS:
%   aep_pos        - [x; y] location of the AEP
%   mu,beta,alpha,delta - system/sail parameters
%
% OUTPUTS:
%   is_stable   - true if all 4 eigenvalues have Re(lambda) <= tol (linearly stable)
%   eigenvalues - 4x1 vector of eigenvalues of the system matrix A
%   J           - 2x2 "potential Hessian" [Uxx Uxy; Uxy Uyy] for reference

    x0 = aep_pos(1);
    y0 = aep_pos(2);

    % Total acceleration field (gravity + sail), planar (z=0), as a
    % function handle of [x,y] for finite-difference Hessian computation.
    accel = @(x,y) totalAccelPlanar(x, y, mu, beta, alpha, delta);

    h = 1e-6; % finite-difference step (nondimensional units)

    % Central differences for second partial derivatives of the potential.
    % Note: Omega_x(x,y) plays the role of dPhi/dx, so
    %   Uxx = d(Omega_x)/dx,  Uxy = d(Omega_x)/dy = d(Omega_y)/dx,  Uyy = d(Omega_y)/dy

    [ax_xp, ay_xp] = accel(x0+h, y0);
    [ax_xm, ay_xm] = accel(x0-h, y0);
    [ax_yp, ay_yp] = accel(x0, y0+h);
    [ax_ym, ay_ym] = accel(x0, y0-h);

    Uxx = (ax_xp - ax_xm) / (2*h);
    Uyy = (ay_yp - ay_ym) / (2*h);
    Uxy = (ax_yp - ax_ym) / (2*h);   % = d(ax)/dy
    Uyx = (ay_xp - ay_xm) / (2*h);   % = d(ay)/dx  (should ~= Uxy if field is conservative-like)

    Uxy_sym = 0.5*(Uxy + Uyx); % symmetrize to suppress FD asymmetry noise

    J = [Uxx, Uxy_sym; Uxy_sym, Uyy];

    % 4x4 linearized system matrix (planar CR3BP variational form)
    A = [0,    0,     1,  0;
         0,    0,     0,  1;
         Uxx,  Uxy_sym, 0,  2;
         Uxy_sym, Uyy,  -2, 0];

    eigenvalues = eig(A);

    tol = 1e-8;
    is_stable = all(real(eigenvalues) <= tol);
end

function [ax, ay] = totalAccelPlanar(x, y, mu, beta, alpha, delta)
% Helper: total planar acceleration (Omega_x + a_sail_x, Omega_y + a_sail_y)
    z = 0;
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

    Omega_x = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    Omega_y = y - (1-mu)*y/r1^3 - mu*y/r2^3;

    a_sail = solarSailAccel(x, y, z, mu, beta, alpha, delta);

    ax = Omega_x + a_sail(1);
    ay = Omega_y + a_sail(2);
end
