function [X,Y,Z,A1,A2,SigE] = make_XY_once(n,p,q,m,scenario)
% scenario:
%   'indep_all'                 : Setting S1 in Section R.3 of the supplementary material  
%   'indep_vecs_corr_within'    : Setting S2 in Section R.3 of the supplementary material  
%   'corr_between_indep_within' : Setting S3 in Section R.3 of the supplementary material 
%   'corr_all'                  : Setting S4 in Section R.3 of the supplementary material 

% X = Z*A1' + E1,  Y = Z*A2' + E2
% Z ~ N(0, I_m);  [E1,E2] from gen_gauss_XY

    % generate noise (E_1,E_2) according to Settings S1-S4 (depends on "scenario")
    % SigE is the covariance matrix of (E1,E2) 
    [E1,E2,SigE] = gen_gauss_XY(n,p,q,scenario);    

    % generate factors and loadings in Equation (R.1) of the supplementary material 
    Z  = randn(n,m);            % n-by-m matrix 
    A1 = rand(p,m);              
    A2 = rand(q,m);             

    % construct X = Z*A1' + E1,  Y = Z*A2' + E2
    X = Z*A1.' + E1;            % n-by-p matirx
    Y = Z*A2.' + E2;            % n-by-q matrix
end
