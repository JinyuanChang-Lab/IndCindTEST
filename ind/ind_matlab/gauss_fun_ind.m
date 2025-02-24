function sA = gauss_fun_ind(x,y, n,p,q, kap1, alpha)
%% generate data
sA = NaN(1, length(kap1)*3);

for s= 1:length(kap1)
    kap =kap1(s);
    x_norm = data_transform(x,kap);
    y_norm = data_transform(y,kap);
    
    %% critival value
    N=5000;
    M=randsrc(N, n);
    M_norm = randn(N, n);
    M_mm =reshape(rmemmer(N*n), N,n);
    tt=100;
    ind_my = p/tt;
    dat_max= zeros(N, ind_my);
    dat_max_n= zeros(N, ind_my);
    dat_max_mm= zeros(N, ind_my);

    for i= 1 : ind_my
       ind1 = (i-1)*tt+1;
       ind2 = i*tt;
       dat = reshape(repmat(x_norm(:, ind1:ind2), q, 1),n,tt*q)  .* repmat(y_norm, 1, tt);
       dat = dat - repmat(mean(dat, 1), n, 1);
      %%%%%radmecher
     %dat_new = abs(M * dat ./ sqrt(n));
       dat_max(:,i) = max(abs(M * dat ./ sqrt(n)), [], 2);
     %%%%%gaussian
     %dat_new_n = abs(M_norm * dat ./ sqrt(n));
       dat_max_n(:,i) = max(abs(M_norm * dat ./ sqrt(n)), [], 2);
     %%%%%memmer
     %dat_new_mm = abs(M_mm * dat ./ sqrt(n));
       dat_max_mm(:,i) = max(abs(M_mm * dat ./ sqrt(n)), [], 2);
    end
%%%%%radmecher
    dat_max1 = sort(max(dat_max, [], 2) , 'descend');
    cv_red = dat_max1(floor(N*alpha),1);
%%%%%gaussian
    dat_max_n1 = sort(max(dat_max_n, [], 2) , 'descend');
    cv_norm = dat_max_n1(floor(N*alpha),1);
%%%%%memmer
    dat_max_mm1 = sort(max(dat_max_mm, [], 2) , 'descend');
    cv_mm = dat_max_mm1(floor(N*alpha),1);

    %%test_sta
    test_1 = sqrt(n) .* max(max(  abs(  x_norm.' * y_norm ./n ) ));

    sA(1,(s-1)*3+1) = (test_1 > cv_red);
    sA(1,(s-1)*3+2) = (test_1 > cv_norm);
    sA(1,(s-1)*3+3) = (test_1 > cv_mm);

   % disp(kap);
end


end

