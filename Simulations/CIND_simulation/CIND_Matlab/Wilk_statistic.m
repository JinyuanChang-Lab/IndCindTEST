function W = Wilk_statistic(X, p1, p2)
    %% the classical likelihood tests (LRT) in Jiang and Yang (2013),
    [n, p] = size(X);

    %   S = 1/n * sum((xi - x_bar) * (xi - x_bar)')
    x_bar = mean(X, 1);   
    S = zeros(p, p);   
    for i = 1:n
        xi = X(i, :);
        S = S + (xi - x_bar)' * (xi - x_bar);   
    end
    S = S / n;  % S = 1/n * sum

    % block S
    S11 = S(1:p1, 1:p1);      
    S12 = S(1:p1, p1+1:end);  
    S21 = S(p1+1:end, 1:p1);  
    S22 = S(p1+1:end, p1+1:end);  

     
    % % Convert to determinants of scatter matrices A = (n-1)S, A11 = (n-1)S11, A22 = (n-1)S22
    % det_A corresponds to |A|, det_A11 to |A11|, det_A22 to |A22|
    
    det_A = (n-1)^p * det(S); 
    term_1 = det_A ;

    %  
    det_A11 = (n-1)^(p1) * det(S11);
    det_A22 = (n-1)^(p2) * det(S22);
    term_2 = (det_A11 * det_A22) ;

    %  % Wilks' Lambda based on block determinants: ¦«_n = |A| / (|A11| |A22|)
    W_n =  term_1 / term_2 ;
    
    rn = (-log(1-p/(n-1)))^(1/2);
    rn1 = (-log(1-p1/(n-1)))^(1/2);
    rn2 = (-log(1-p2/(n-1)))^(1/2);
    %  
    mun = -rn^2 * (p - (n-1) + 1/2) + rn1^2 * (p1 - (n-1) + 1/2) + rn2^2 * (p2 - (n-1) + 1/2);
    %  
    sigma_n2 = 2 * rn^2 - 2 * (rn1^2 + rn2^2);
 
    W = (log(W_n) - mun) /sigma_n2^(1/2);
    return;
       
end
