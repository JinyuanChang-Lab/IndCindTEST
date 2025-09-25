function [X, Y] = generate_gauss_data(n, p, q)
    % Generate Gaussian data from Setting S5 in Section R.3 of the supplementary material. 
    % Inputs:
    %   n : sample size
    %   p : dimension of X
    %   q : dimension of Y
    
    % Build covariance matrix Sigma = 0.75 * eye(p+q) + 0.25 * ones(p+q)
    Sigma = 0.75 * eye(p + q) + 0.25 * ones(p + q);
    
    % Draw standard normal samples
    Z = randn(n, p + q);  % each row ~ N(0, I_{p+q})
    
    % Check whether Sigma is positive definite.
    [L, flag] = chol(Sigma, 'lower');
    if flag ~= 0
        error('matrix degenerate');
    end
    
    % Generate Gaussian data with covariance matrix Sigma
    W = Z * L.';   
    
    % Split into X (first p columns) and Y (last q columns)
    X = W(:, 1:p);
    Y = W(:, p+1:end);
   
end

