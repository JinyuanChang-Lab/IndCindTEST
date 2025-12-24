function [W_w, W_lh, W_bnp] = TB_statistics(X, p1, p2)
    %% three different tests in Bodnar et al. (2019),
    % X is an n x p matrix
    % p1 and p2 are the partition sizes, with p1 + p2 = p
    % k is the number of partitions (k = 2 in this case)

    [n, p] = size(X);
    
    % Partition X into two blocks based on p1 and p2
    X1 = X(:, 1:p1);   % First p1 columns
    X2 = X(:, p1+1:end);  % Last p2 columns

    % Calculate the sample covariance matrix S
    S = (X' * X) / n; % Sample covariance matrix

    % Partition the sample covariance matrix
    S11 = n* S(1:p1, 1:p1);
    S12 = n* S(1:p1, p1+1:end);
    S21 = n* S(p1+1:end, 1:p1);
    S22 = n* S(p1+1:end, p1+1:end);

    % Compute the inverses of S11 and S22
    S11_inv = inv(S11);

    % Compute W and T
    a = S21 * S11_inv * S12;
    W = a / p1;
    T = (S22 - a) /(n-p1);

    % Compute the eigenvalues of W and T
    eig_WT_inv = eig(W * inv(T));  
    
    
    eig_WT_inv = sort(eig_WT_inv, 'descend');   
    
    %  Wilks' Lambda  
    TW=0;  TLH=0; TBNP =0;
    for i = 1:(p-p1)
        TW = TW +  log(1 + eig_WT_inv(i)) ;  % Wilks' Lambda
        TLH = TLH +  eig_WT_inv(i);  % Lawley-Hotelling
        TBNP = TBNP + eig_WT_inv(i)/(1+eig_WT_inv(i));  % Bartlett-Nanda-Pillai
    end
    
    %   gamma_1n and gamma_2n
    gamma_1n = (p - p1) / p1;  %   gamma_1n
    gamma_2n = (p - p1) / (n - p1);  %   gamma_2n
    
    %  h_n
    h_n = sqrt(gamma_1n + gamma_2n - gamma_1n * gamma_2n);  %  h_n
    
    [w_n, d_n] = solve_wd(gamma_2n, h_n);
    
    %   mu_w
    mu_w = 0.5 * log(( (w_n^2 - d_n^2) * h_n^2 ) / ( (w_n * h_n - gamma_2n * d_n)^2 ));
    
    %   sigma_w^2
    sigma_w  = 2 * log(w_n^2 / (w_n^2 - d_n^2));
    
    %   mu_LH and sigma_LH^2
    mu_LH = gamma_2n / (1 - gamma_2n)^2;  %   mu_LH
    sigma_LH = 2 * h_n^2 / (1 - gamma_2n)^4;  %   sigma_LH^2
    
    %   mu_BNP and sigma_BNP^2
    mu_BNP = (1 - gamma_2n)^2 * w_n^2 * (d_n^2 - gamma_2n) / (w_n^2 - d_n^2)^2;  %   mu_BNP
    sigma_BNP = 2 * d_n^2 * (1 - gamma_2n)^4 * (w_n^2 * (w_n + d_n) + d_n^3 * (w_n^2 - 1)) / ...
                (w_n^2 * (1 + d_n) * (w_n^2 - d_n^2)^4);  %   sigma_BNP^2
            
    %   s_w 
    b = -log((1 - gamma_2n)^2) - (1 - gamma_2n) / gamma_2n * log(w_n) + ...
              (gamma_1n + gamma_2n) / (gamma_1n * gamma_2n) * log(w_n - d_n * gamma_2n / h_n);
    if gamma_1n > 1
        s_w = b - (1-gamma_1n) / gamma_1n * log(w_n -d_n / h_n);
    elseif gamma_1n == 1
        s_w = b;
    else
        s_w = b - (1-gamma_1n) / gamma_1n * log(w_n -d_n * h_n);
    end
    
    %   s_LH and s_BNP
    s_LH = 1 / (1 - gamma_2n);  %   s_LH
    s_BNP = 1 / (w_n^2 - gamma_2n);  %   s_BNP     
    
    W_w = (TW - (p-p1) * s_w -mu_w ) / sigma_w^(1/2);
    W_lh = (TLH - (p-p1) * s_LH -mu_LH ) / sigma_LH^(1/2);
    W_bnp = (TBNP - (p-p1) * s_BNP -mu_BNP ) / sigma_BNP^(1/2);
     
end

%% help function

function [w_n, d_n] = solve_wd(gamma_2n, h_n)
    %% INPUT:
    %   gamma_2n : given scalar gamma_2n
    %   h_n      : given scalar h_n
    %
    % OUTPUT:
    %   w_n, d_n : solution pair satisfying w_n^2 + (h_n^2)/(w_n^2) = constant_term
    %              with constant_term = (1 - gamma_2n)^2 + 1 + h_n^2
    rng(1);
    % Compute the constant term on the right-hand side
    constant_term = (1 - gamma_2n)^2 + 1 + h_n^2;
    
    % Coefficients of the quartic in w_n arranged as a quadratic in w_n^2:
    %   w_n^2 + h_n^2 / w_n^2 = constant_term
    % => w_n^4 - constant_term * w_n^2 + h_n^2 = 0
    
    a = 1;                 % coefficient of w_n^4
    b = -constant_term;    % coefficient of w_n^2
    c = h_n^2;             % constant term
    
    % Solve the quadratic: a * (w_n^2)^2 + b * (w_n^2) + c = 0
    discriminant = b^2 - 4 * a * c;  % discriminant
    
    if discriminant >= 0
        % Compute the two candidates for w_n^2
        w_n_squared_1 = (-b + sqrt(discriminant)) / (2 * a);
        w_n_squared_2 = (-b - sqrt(discriminant)) / (2 * a);
        
        % Choose the valid (nonnegative) solution for w_n^2
        w_n_squared = max(w_n_squared_1, w_n_squared_2);
        
        % Recover w_n and d_n
        w_n = sqrt(w_n_squared);   % w_n >= 0
        d_n = h_n / w_n;           % corresponding d_n
    else
        error('No real solution exists!');
    end
end
