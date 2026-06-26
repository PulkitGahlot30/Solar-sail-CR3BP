function dataset = generateDataset(varargin)
% GENERATEDATASET  Sweeps the solar sail CR3BP parameter space
% (mu, beta, alpha, delta), locates the AEP near L4 (and L5 by symmetry),
% performs linear stability analysis at each point, and assembles a
% table suitable for training a PINN / ML stability classifier.
%
% Optional name-value args:
%   'MuRange'     - [mu_min, mu_max, n_mu]      default [0.001, 0.5, 12]
%   'BetaRange'   - [beta_min, beta_max, n_beta] default [0.0, 0.30, 16]
%   'AlphaRange'  - [alpha_min, alpha_max, n_alpha] (deg) default [0, 60, 13]
%   'DeltaRange'  - [delta_min, delta_max, n_delta] (deg) default [0, 0, 1]
%   'OutFile'     - output CSV path, default 'data/stability_dataset.csv'
%
% OUTPUT:
%   dataset - MATLAB table with columns:
%     mu, beta, alpha_deg, delta_deg, x_aep, y_aep, converged,
%     lambda1_re, lambda1_im, lambda2_re, lambda2_im,
%     lambda3_re, lambda3_im, lambda4_re, lambda4_im,
%     max_real_eig, is_stable

    p = inputParser;
    addParameter(p, 'MuRange', [0.001, 0.45, 12]);
    addParameter(p, 'BetaRange', [0.0, 0.30, 16]);
    addParameter(p, 'AlphaRange', [0, 60, 13]);   % degrees
    addParameter(p, 'DeltaRange', [0, 0, 1]);     % degrees (planar case: delta=0 default)
    addParameter(p, 'OutFile', 'data/stability_dataset.csv');
    parse(p, varargin{:});
    R = p.Results;

    mu_vals    = linspace(R.MuRange(1), R.MuRange(2), R.MuRange(3));
    beta_vals  = linspace(R.BetaRange(1), R.BetaRange(2), R.BetaRange(3));
    alpha_vals = linspace(R.AlphaRange(1), R.AlphaRange(2), R.AlphaRange(3));
    delta_vals = linspace(R.DeltaRange(1), R.DeltaRange(2), R.DeltaRange(3));

    n_total = numel(mu_vals)*numel(beta_vals)*numel(alpha_vals)*numel(delta_vals);
    fprintf('Total parameter combinations to evaluate: %d\n', n_total);

    % Preallocate storage
    rows = cell(n_total, 1);
    idx = 0;
    t_start = tic;

    for im = 1:numel(mu_vals)
        mu = mu_vals(im);
        % Classical L4 as the continuation seed (sail perturbs it smoothly
        % for small/moderate beta)
        L4_classical = [0.5 - mu, sqrt(3)/2];

        for ib = 1:numel(beta_vals)
            beta = beta_vals(ib);
            % Use previous beta's converged AEP as warm-start when possible
            guess = L4_classical;

            for ia = 1:numel(alpha_vals)
                alpha = deg2rad(alpha_vals(ia));

                for id = 1:numel(delta_vals)
                    delta = deg2rad(delta_vals(id));
                    idx = idx + 1;

                    [aep_pos, converged] = findAEP(mu, beta, alpha, delta, guess);

                    if converged
                        guess = aep_pos; % warm-start next iteration
                        [is_stable, eigvals, ~] = stabilityAnalysis(aep_pos, mu, beta, alpha, delta);
                        max_re = max(real(eigvals));
                        lam_re = real(eigvals);
                        lam_im = imag(eigvals);
                    else
                        is_stable = false;
                        eigvals = NaN(4,1);
                        max_re = NaN;
                        lam_re = NaN(4,1);
                        lam_im = NaN(4,1);
                        aep_pos = [NaN; NaN];
                    end

                    rows{idx} = {mu, beta, rad2deg(alpha), rad2deg(delta), ...
                        aep_pos(1), aep_pos(2), converged, ...
                        lam_re(1), lam_im(1), lam_re(2), lam_im(2), ...
                        lam_re(3), lam_im(3), lam_re(4), lam_im(4), ...
                        max_re, is_stable};
                end
            end
        end
        fprintf('mu = %.4f done (%d/%d), elapsed %.1f s\n', mu, im, numel(mu_vals), toc(t_start));
    end

    varnames = {'mu','beta','alpha_deg','delta_deg','x_aep','y_aep','converged', ...
        'lambda1_re','lambda1_im','lambda2_re','lambda2_im', ...
        'lambda3_re','lambda3_im','lambda4_re','lambda4_im', ...
        'max_real_eig','is_stable'};

    dataset = cell2table(vertcat(rows{:}), 'VariableNames', varnames);

    % Ensure output directory exists
    [outdir, ~, ~] = fileparts(R.OutFile);
    if ~isempty(outdir) && ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    writetable(dataset, R.OutFile);
    fprintf('Dataset written to %s (%d rows)\n', R.OutFile, height(dataset));
end
