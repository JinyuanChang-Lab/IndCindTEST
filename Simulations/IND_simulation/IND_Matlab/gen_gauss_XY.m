function [X, Y, Sigma] = gen_gauss_XY(n, p, q, mode)
% Generate Gaussian samples (X, Y) under four scenarios:
% mode:
%   'indep_all'                 : X and Y are independent; components within each are independent  (S1 in Section R.3 of the supplementary material) 
%   'indep_vecs_corr_within'    : X and Y are independent; within-vector components are correlated (S2 in Section R.3 of the supplementary material) 
%   'corr_between_indep_within' : X and Y are correlated; within-vector components are independent (S3 in Section R.3 of the supplementary material)
%   'corr_all'                  : X and Y are correlated; within-vector components are correlated  (S4 in Section R.3 of the supplementary material)
%

    % set the correlation parameters
    opts = struct; 
    rhoX   = getdef(opts, 'rhoX',  0.5);
    rhoY   = getdef(opts, 'rhoY',  0.5);
    rhoXY  = getdef(opts, 'rhoXY', 0.5);
    within = getdef(opts, 'within','ar1');

    % --- Within-block covariance (correlation matrix with unit variances) ---
    SigmaX = eye(p);
    SigmaY = eye(q);
    if any(strcmpi(mode, {'indep_vecs_corr_within','corr_all'}))
        SigmaX = corr_mat(p, rhoX, within);
        SigmaY = corr_mat(q, rhoY, within);
    end

    % --- Cross-block covariance C) ---
    C = zeros(p, q);
    if any(strcmpi(mode, {'corr_between_indep_within','corr_all'}))
        r = min(p, q);
        C(1:r, 1:r) = rhoXY * eye(r);  % pair X1~Y1, ..., Xr~Yr with correlation coefficient rhoXY
    end

    % --- The covariance matrix of (X,Y) ---
    Sigma = [SigmaX, C; C', SigmaY];

    % --- If Sigma is not positive definite, shrink the cross-block correlations ---
    [L, flag] = chol(Sigma, 'lower');
    if flag ~= 0
        scale = 0.98;
        while flag ~= 0
            C = C * scale;
            Sigma = [SigmaX, C; C', SigmaY];
            [L, flag] = chol(Sigma, 'lower');
            scale = scale * 0.98;
            if scale < 1e-6
                error('matrix is degenerate');
            end
        end
    end
    
    % --- generate data ---
    Z = randn(n, p+q);     % each row ~ N(0, I_{p+q})
    W = Z * L.';           % impose covariance Sigma on rows
    X = W(:, 1:p);
    Y = W(:, p+1:end);
end

%  within-block correlation matrix (unit variances)
function S = corr_mat(d, rho, type)
    if d == 0, S = []; return; end
    switch lower(type)
        case 'ar1'   % AR(1): S_ij = rho^{|i-j|}, |rho|<1 
            idx = (0:d-1)';
            S = rho .^ abs(idx - idx');
        case 'cs'    % Compound symmetry/equicorrelation: diag=1, off-diag=rho 
            if ~(rho > -1/(d-1) && rho < 1)
                error('rho need (-1/(d-1), 1).');
            end
            S = (1-rho)*eye(d) + rho*ones(d);
        otherwise
            error('unknown');
    end
end

% get struct field with default
function v = getdef(s, f, def)
    if isfield(s, f), v = s.(f); else, v = def; end
end
