function [X,Y,Z,A1,A2] = make_XY_once_dence(n,p,q,m)
% Setting S5 in Section R.3 of the supplementary material
% X = Z*A1' + E1,  Y = Z*A2' + E2
%  Z ~ N(0, I_m);  [E1,E2] from gen_gauss_XY

    [E1,E2] = generate_gauss_data(n, p, q);
    
    % generate factors and loadings in Equation (R.1) of the supplementary material 
    Z  = randn(n,m);            % n-by-m matrix 
    A1 = rand(p,m);             
    A2 = rand(q,m);   

    % construct X = Z*A1' + E1,  Y = Z*A2' + E2  
    X = Z*A1.' + E1;            % n-by-p matrix
    Y = Z*A2.' + E2;            % n-by-q matrix
end