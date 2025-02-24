%% ex1
clear
parpool(120);
ex=2;


n = 50;
p=100; q=100;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;

dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);


n = 100;
p=100; q=100;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;

dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);




n = 50;
p=400; q=400;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;
dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);


n = 100;
p=400; q=400;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;

dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

n = 50;
p=1600; q=1600;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;

dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);


n = 100;
p=1600; q=1600;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)
alpha=0.05;
%kap1=[5/6,6/7,7/8,8/9,9/10];
%%%%%number
n_my = 1000;

dep=0;
res_my = NaN(n_my,  3);
%title = ['myresex',num2str(ex), '_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
tic;
parfor i=1 : n_my
    %tic
    rng(i*100+123);
    [x,y] = data1_ex1(n,p,q,dep);
    res1 = gauss_fun_ind(x,y,n,p,q, alpha);
    res_my(i,:)= res1;
    %toc
    %writetable(table(res_my),title);    
    fprintf('the %d time \n', i)
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);





delete(gcp('nocreate'));