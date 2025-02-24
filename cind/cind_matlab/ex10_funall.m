%% ex10
clear
addpath('/home/yuedu/glmnet_matlab')
parpool(120);
alpha=0.05;

%%%%%number
n_my = 1000;
ex=10;
funcName = sprintf('data1_ex%d', ex);

%% 100-100;
p=100; q=100; m=5;
n=100;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    %[x,y,z] = data1_ex6(n,p,q, m,dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

%% 100-200;
p=100; q=100; m=5;
n=200;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);


%% 400-100;
p=400; q=400; m=5;
n=100;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

%% 400-200;
p=400; q=400; m=5;
n=200;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);


%% 1600-100;
p=1600; q=1600; m=5;
n=100;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

%% 1600-200;
p=1600; q=1600; m=5;
n=200;
fprintf('-----------  ex %d ----------- \n', ex)
fprintf('-----------  n %d ----------- \n', n)
fprintf('-----------  p %d ----------- \n', p)

dep =0;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);
   
 
dep =p/10;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

dep =p/5;
fprintf('-----------  dep %d ----------- \n', dep)
res_my = NaN(n_my,  3);
tic;
for i=1 : n_my
    rng(i*100);
    [x,y,z] =feval(funcName, n, p, q, m, dep);
    res1= gauss_fun(x,y,z,n,p,q,m, alpha);
    res_my(i,:)= res1;
    if mod(i, 200) == 0
         fprintf('the %d time \n', i);
    end
end
toc
result = mean(res_my, 1);
title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
writetable(table(result),title1);

delete(gcp('nocreate'));