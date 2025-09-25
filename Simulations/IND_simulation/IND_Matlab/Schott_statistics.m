function Schott = Schott_statistics(X, p1, p2)
    %% the test statistics of  (Schott) in Bao et al. (2017) 
     
    [n, p] = size(X);
    X1 = (X(:, 1:p1))';
    X2 = (X(:,(p1+1):p))';
    
    
    C21 = (X2*X2')^(-1/2) * (X2*X1') * (X1*X1')^(-1) * (X1*X2') * (X2*X2')^(-1/2);
    %C22 = (X2*X2')^(-1/2) * (X2*X2') * (X2*X2')^(-1) * (X2*X2') * (X2*X2')^(-1/2);
    
    %sB  = (1/2) * (trace(C11) + trace(C12) +trace(C21) + trace(C22)) - p/2;
    
    sB  = trace(C21); 
    % Centering and scaling constants (large-sample approximations)
    an =  (p1*p2) /(n-1);
    bn = 2 * p1 *p2 * (n-1-p1)*(n-1-p2) /(n-1)^4;
    
    Schott = (sB -an) / (bn)^(1/2);
    
    
    
    return;
end
