function [sA,sB] = gauss_data_ind(x,y)
%% generate data
sA = NaN(1, 3);
sB = NaN(1, 3);
[n,p] = size(x);   q = size(y,2); 

x_norm = data_transform_new(x);
y_norm = data_transform_new(y);
    
  
%% critival value
N=5000;
M=randsrc(N, n);
M_norm = randn(N, n);
M_mm =reshape(rmemmer(N*n), N,n);
dat = reshape(repmat(x_norm , q, 1),n,p*q).* repmat(y_norm , 1, p);
dat = dat - repmat(mean(dat, 1), n, 1);
    %%%%%radmecher
    %dat_new = abs(M * dat ./ sqrt(n));
dat_max = sort(max(abs(M * dat ./ sqrt(n)), [], 2) , 'descend');
%cv_red = dat_max(N*alpha,1);
    %%%%%gaussian
    %dat_new_n = abs(M_norm * dat ./ sqrt(n));
dat_max_n = sort(max(abs(M_norm * dat ./ sqrt(n)), [], 2) , 'descend');
%cv_norm = dat_max_n(N*alpha,1);
    %%%%%memmer
    %dat_new_mm = abs(M_mm * dat ./ sqrt(n));
dat_max_mm = sort(max(abs(M_mm * dat ./ sqrt(n)), [], 2) , 'descend');
%cv_mm = dat_max_mm(N*alpha,1);
    
   
%% 
    %%%%%%
test_1 = sqrt(n) * max(max(abs(x_norm .' * y_norm ./n)));
sA(1,1) = mean(dat_max > test_1);
sB(1,1) = test_1;
    %%%%%%
sA(1,2) = mean(dat_max_n > test_1);
sB(1,2) = test_1;   
    %%%%%%
sA(1,3) = mean(dat_max_mm > test_1);
sB(1,3) = test_1;


end


