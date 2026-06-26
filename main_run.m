%% main_run.m
% Master script for the Solar Sail CR3BP Triangular AEP Stability pipeline.
% Run sections in order (Ctrl+Enter per cell in MATLAB, or run whole file).

clear; clc; close all;

%% --- 0. Sanity check: classical (beta=0) case must recover textbook L4 result ---
% For beta = 0 (no sail), the AEP must coincide with the classical L4
% point, and stability must match the well-known result:
%   stable for mu < mu_c = (1 - sqrt(23/27))/2 ~ 0.0385

mu_test = 0.01;  % well inside the stable regime
beta0 = 0; alpha0 = 0; delta0 = 0;
L4_guess = [0.5 - mu_test, sqrt(3)/2];

[aep0, conv0] = findAEP(mu_test, beta0, alpha0, delta0, L4_guess);
fprintf('--- Sanity check (mu=%.4f, beta=0) ---\n', mu_test);
fprintf('AEP found at: (%.8f, %.8f)\n', aep0(1), aep0(2));
fprintf('Classical L4: (%.8f, %.8f)\n', L4_guess(1), L4_guess(2));
fprintf('Converged: %d\n', conv0);

[stab0, eig0, ~] = stabilityAnalysis(aep0, mu_test, beta0, alpha0, delta0);
fprintf('Is stable: %d\n', stab0);
fprintf('Eigenvalues:\n'); disp(eig0);
mu_c = (1 - sqrt(23/27))/2;
fprintf('Critical mass ratio mu_c = %.6f (mu_test should be stable since mu_test < mu_c: %d)\n\n', mu_c, mu_test < mu_c);

% Quick check just above mu_c -> should be unstable
mu_test2 = 0.045;
L4_guess2 = [0.5 - mu_test2, sqrt(3)/2];
[aep0b, ~] = findAEP(mu_test2, beta0, alpha0, delta0, L4_guess2);
[stab0b, eig0b, ~] = stabilityAnalysis(aep0b, mu_test2, beta0, alpha0, delta0);
fprintf('--- Sanity check (mu=%.4f > mu_c, beta=0) ---\n', mu_test2);
fprintf('Is stable: %d  (expected: 0, unstable)\n\n', stab0b);

%% --- 1. Small smoke-test sweep (fast, for debugging) ---
fprintf('Running SMALL smoke-test sweep...\n');
small_data = generateDataset( ...
    'MuRange',   [0.01, 0.04, 4], ...
    'BetaRange', [0.0, 0.1, 4], ...
    'AlphaRange',[0, 30, 4], ...
    'DeltaRange',[0, 0, 1], ...
    'OutFile',   'data/smoke_test.csv');

disp(small_data(1:min(10,height(small_data)), :));

%% --- 2. Visualize smoke test ---
plotStabilityMap(small_data, 0.02);
saveas(gcf, 'data/smoke_test_stability_map.png');

%% --- 3. FULL parameter sweep for PINN training dataset ---
% WARNING: depending on grid resolution this can take a while since each
% point calls fsolve + finite-difference Hessian. Tune resolution as needed.
fprintf('Running FULL sweep for PINN training dataset...\n');
full_data = generateDataset( ...
    'MuRange',   [0.001, 0.45, 30], ...
    'BetaRange', [0.0, 0.30, 30], ...
    'AlphaRange',[0, 60, 25], ...
    'DeltaRange',[0, 0, 1], ...
    'OutFile',   'data/stability_dataset.csv');

fprintf('Full dataset size: %d rows\n', height(full_data));
fprintf('Stable fraction: %.2f%%\n', 100*mean(full_data.is_stable(full_data.converged)));

%% --- 4. Visualize full stability maps for a few representative mu values ---
mu_show = [0.01, 0.02, 0.03, 0.0385];
for k = 1:numel(mu_show)
    plotStabilityMap(full_data, mu_show(k));
    saveas(gcf, sprintf('data/stability_map_mu_%0.4f.png', mu_show(k)));
end

fprintf('Pipeline complete. Dataset and figures saved in ./data/\n');
