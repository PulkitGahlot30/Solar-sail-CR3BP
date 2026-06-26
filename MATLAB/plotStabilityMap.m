function plotStabilityMap(dataset, mu_fixed)
% PLOTSTABILITYMAP  Visualizes the (beta, alpha) stability map for a
% fixed mass ratio mu, using the dataset produced by generateDataset.m.
%
% INPUTS:
%   dataset   - table from generateDataset.m (or readtable on the CSV)
%   mu_fixed  - the mu value to slice on (nearest match used)

    mu_vals = unique(dataset.mu);
    [~, idx] = min(abs(mu_vals - mu_fixed));
    mu_use = mu_vals(idx);

    sub = dataset(abs(dataset.mu - mu_use) < 1e-9, :);

    figure('Color', 'w', 'Position', [100 100 700 550]);
    hold on;

    stable_rows   = sub.converged & sub.is_stable;
    unstable_rows = sub.converged & ~sub.is_stable;
    nonconv_rows  = ~sub.converged;

    scatter(sub.beta(stable_rows), sub.alpha_deg(stable_rows), 45, 'g', 'filled', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Stable AEP');
    scatter(sub.beta(unstable_rows), sub.alpha_deg(unstable_rows), 45, 'r', 'filled', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Unstable AEP');
    scatter(sub.beta(nonconv_rows), sub.alpha_deg(nonconv_rows), 45, [0.6 0.6 0.6], 'x', ...
        'DisplayName', 'No convergence');

    xlabel('Sail lightness number, \beta');
    ylabel('Pitch angle, \alpha (deg)');
    title(sprintf('Stability map of AEP near L_4 (\\mu = %.4f)', mu_use));
    legend('Location', 'best');
    grid on;
    box on;
    hold off;
end
