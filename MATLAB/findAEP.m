function [aep_pos, converged] = findAEP(mu, beta, alpha, delta, initial_guess)
% FINDAEP  Locates an Artificial Equilibrium Point (AEP) near the
% triangular libration point using fsolve, starting from a supplied
% initial guess (typically the classical L4 or L5 position).
%
% INPUTS:
%   mu             - CR3BP mass ratio
%   beta           - sail lightness number
%   alpha, delta   - sail pitch/clock angles [rad]
%   initial_guess  - [x0; y0], e.g. classical L4 = [0.5-mu, sqrt(3)/2]
%
% OUTPUTS:
%   aep_pos     - [x; y] location of the AEP (NaN if not converged)
%   converged   - logical flag

    opts = optimoptions('fsolve', 'Display', 'off', ...
        'TolFun', 1e-13, 'TolX', 1e-13, 'MaxIterations', 500);

    try
        [sol, fval, exitflag] = fsolve(@(v) aepResidual(v, mu, beta, alpha, delta), ...
            initial_guess, opts);
        converged = (exitflag > 0) && (norm(fval) < 1e-9);
        if converged
            aep_pos = sol;
        else
            aep_pos = [NaN; NaN];
        end
    catch
        aep_pos = [NaN; NaN];
        converged = false;
    end
end
