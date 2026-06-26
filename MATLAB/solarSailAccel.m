function a_sail = solarSailAccel(x, y, z, mu, beta, alpha, delta)
% SOLARSAILACCEL  Solar Radiation Pressure acceleration vector for a flat
% solar sail in the Circular Restricted Three-Body Problem (CR3BP),
% following the standard McInnes (1999) sail force model.
%
% INPUTS:
%   x,y,z   - position of the sail in the synodic (rotating) frame
%             (nondimensional, primaries at (-mu,0,0) and (1-mu,0,0))
%   mu      - mass ratio of the CR3BP system (m2/(m1+m2))
%   beta    - sail lightness number (0 <= beta <= 1)
%   alpha   - pitch angle [rad]: angle between sail normal n_hat and
%             the Sun-sail line r_s_hat
%   delta   - clock angle [rad]: angle of n_hat about r_s_hat,
%             measured in the (p,q) plane orthogonal to r_s_hat
%
% OUTPUT:
%   a_sail  - [3x1] sail acceleration vector in synodic frame
%
% NOTE: The larger primary (mass 1-mu, the "Sun") is assumed to be
% located at (-mu, 0, 0). The sail responds to radiation from this body.

    % Sun-to-sail vector and distance
    rs_vec = [x + mu; y; z];      % vector from Sun (at -mu,0,0) to sail
    rs = norm(rs_vec);
    if rs < 1e-12
        error('solarSailAccel: sail position coincides with primary (rs ~ 0).');
    end
    rs_hat = rs_vec / rs;

    % Build orthonormal frame {rs_hat, p_hat, q_hat} following McInnes:
    % p_hat lies in the orbital (x-y) plane, orthogonal to rs_hat
    % q_hat completes the right-handed triad
    z_axis = [0; 0; 1];
    p_hat = cross(z_axis, rs_hat);
    if norm(p_hat) < 1e-9
        % rs_hat nearly parallel to z-axis; pick an alternate reference
        p_hat = cross([0;1;0], rs_hat);
    end
    p_hat = p_hat / norm(p_hat);
    q_hat = cross(rs_hat, p_hat);
    q_hat = q_hat / norm(q_hat);

    % Sail normal vector n_hat parameterized by pitch (alpha) and clock (delta)
    % n_hat = cos(alpha)*rs_hat + sin(alpha)*cos(delta)*p_hat + sin(alpha)*sin(delta)*q_hat
    n_hat = cos(alpha)*rs_hat + sin(alpha)*cos(delta)*p_hat + sin(alpha)*sin(delta)*q_hat;
    n_hat = n_hat / norm(n_hat);

    % Sail must face the Sun: enforce n_hat . rs_hat >= 0
    cosA = dot(n_hat, rs_hat);
    if cosA < 0
        n_hat = -n_hat;
        cosA = -cosA;
    end

    % McInnes flat-sail SRP acceleration magnitude factor
    % a_sail = beta * (1-mu)/rs^2 * (n_hat . rs_hat)^2 * n_hat
    factor = beta * (1 - mu) / rs^2 * cosA^2;
    a_sail = factor * n_hat;
end
