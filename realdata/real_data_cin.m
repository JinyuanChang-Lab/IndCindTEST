clear
addpath('/home/yuedu/matlab_sim/glmnet_matlab')
parpool(45);
%kap1=[6/7,7/8,8/9,9/10];
alpha=0.05;
kap1=[5/6,9/10];

name = ["CS"  "CD" "CSt"  "Eng"  "Fin" "HC"  "Ind"  "IT"  "Mat"  "RE"  "Uti"];
com_my = nchoosek(name,2);
rep_num = size(com_my,1);
res = NaN(rep_num, 4+6*length(kap1));
for i=1:rep_num
    rng(12345);
    x_name = com_my(i,1);
    y_name = com_my(i,2);
    index1 = find(name == x_name);
    index2 = find(name == y_name);
    z_name =name;
    z_name([index1,index2])=[];
    title_x = strcat('02return_',x_name, '.csv');
    x = importdata(title_x).data;
    %xx=csvread(title_x,1,1);
    
    title_y = strcat('02return_',y_name, '.csv');
    y = importdata(title_y).data;
    t = size(z_name, 2);
    title_z = strcat('02return_',z_name(1), '.csv');
    z = importdata(title_z).data;
    for j=2:t
        title_zz = strcat('02return_',z_name(j), '.csv');
        zz = importdata(title_zz).data;
        z =[z,zz];
    end
   [n,p] = size(x);   q = size(y,2);  m = size(z, 2);
   [sA,sB, size1, n_c] = gauss_data(x,y,z, kap1, alpha);
   res(i,1)= n; res(i,2)=p;   res(i,3)=q;  res(i,4)=m;
   res(i, 5:(5+3*length(kap1)-1)) =sA;
   res(i, (5+3*length(kap1)):(4+6*length(kap1))) =sB;
   %a = strcat(x_name,'-',y_name);
   %res(i,1)=a;
end


delete(gcp('nocreate'));