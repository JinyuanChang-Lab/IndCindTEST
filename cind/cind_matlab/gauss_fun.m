function  sA = gauss_fun(x,y,z, n,p,q,m, alpha)
%% generate data
sA = NaN(1,  3);


x_norm = data_transform_new(x);
y_norm = data_transform_new(y);
z_norm = data_transform_new(z);
    
%% lasso regression
lambda_min_u =NaN(1, p);
lambda_min_v =NaN(1, q);
a_hat =NaN(m,p);
b_hat =NaN(m,q);
    % tic;
    % ticBytes(gcp);
parfor i= 1 : p
    fit1 = cvglmnet(z_norm, x_norm(:, i));
    lambda_min_u(1,i) = fit1.lambda_min;
    a_hat(:, i) = fit1.glmnet_fit.beta(:, find(fit1.lambda==fit1.lambda_min));
    fit2 = cvglmnet(z_norm, y_norm(:, i));
    lambda_min_v(1,i) = fit2.lambda_min;
    b_hat(:, i) = fit2.glmnet_fit.beta(:, find(fit2.lambda==fit2.lambda_min));
end
eps = x_norm - z_norm * a_hat;
delta = y_norm - z_norm * b_hat;
    % toc
    % tocBytes(gcp)
    
    %% critival value
N=5000;
M=randsrc(N, n);
M_norm = randn(N, n);
M_mm =reshape(rmemmer(N*n), N,n);
dat = reshape(repmat(eps, q, 1),n,p*q).* repmat(delta, 1, p);
dat = dat - repmat(mean(dat, 1), n, 1);
    %%%%%radmecher
    %dat_new = abs(M * dat ./ sqrt(n));
dat_max = sort(max(abs(M * dat ./ sqrt(n)), [], 2) , 'descend');
cv_red = dat_max(N*alpha,1);
    %%%%%gaussian
    %dat_new_n = abs(M_norm * dat ./ sqrt(n));
dat_max_n = sort(max(abs(M_norm * dat ./ sqrt(n)), [], 2) , 'descend');
cv_norm = dat_max_n(N*alpha,1);
    %%%%%memmer
    %dat_new_mm = abs(M_mm * dat ./ sqrt(n));
dat_max_mm = sort(max(abs(M_mm * dat ./ sqrt(n)), [], 2) , 'descend');
cv_mm = dat_max_mm(N*alpha,1);
    

%%  

test_1 = sqrt(n) * max(max(abs(eps.' * delta./ n)));
sA(1,1) = (test_1 > cv_red);
    %%%%%%
test_2 = sqrt(n) * max(max(abs(eps.' * delta./n)));
sA(1,2) = (test_2 > cv_norm);
    %%%%%%
test_3 = sqrt(n) * max(max(abs(eps.' * delta./n)));
sA(1,3) = (test_3 > cv_mm);



end

