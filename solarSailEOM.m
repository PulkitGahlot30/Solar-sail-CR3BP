function dstate = solarSailEOM(t, state, mu, beta, alpha, delta)
% SOLARSAILEOM  Full 3D equations of motion for the solar sail CR3BP.
% State vector: [x; y; z; vx; vy; vz]  (synodic/rotating frame, nondimensional)
%
% Equations (McInnes 1999 sail force + classical CR3BP gravity/Coriolis):
%   xdd - 2*yd = Omega_x + a_sail_x
%   ydd + 2*xd = Omega_y + a_sail_y
%   zdd        = Omega_z + a_sail_z
%
% where Omega is the augmented potential of the classical CR3BP:
%   Omega = 1/2*(x^2+y^2) + (1-mu)/r1 + mu/r2

    x = state(1); y = state(2); z = state(3);
    vx = state(4); vy = state(5); vz = state(6);

    r1 = sqrt((x+mu)^2 + y^2 + z^2);       % distance from larger primary (Sun)
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);     % distance from smaller primary

    % Partial derivatives of the classical augmented potential Omega
    Omega_x = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    Omega_y = y - (1-mu)*y/r1^3 - mu*y/r2^3;
    Omega_z =   - (1-mu)*z/r1^3 - mu*z/r2^3;

    % Solar sail acceleration
    a_sail = solarSailAccel(x, y, z, mu, beta, alpha, delta);

    xdd = 2*vy + Omega_x + a_sail(1);
    ydd = -2*vx + Omega_y + a_sail(2);
    zdd = Omega_z + a_sail(3);

    dstate = [vx; vy; vz; xdd; ydd; zdd];
end
