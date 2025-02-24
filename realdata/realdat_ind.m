clear
addpath('/home/yuedu/matlab_sim/glmnet_matlab')
parpool(120);
alpha=0.05;

name = ["CS"  "CD" "CSt"  "Eng"  "Fin" "HC"  "Ind"  "IT"  "Mat"  "RE"  "Uti"];
com_my = nchoosek(name,2);
rep_num = size(com_my,1);
res = NaN(rep_num, 7);
%res_string = string(res);

for i=1:rep_num
    rng(12345);
    x_name = com_my(i,1);
    y_name = com_my(i,2);
    index1 = find(name == x_name);
    index2 = find(name == y_name);
    title_x = strcat('02return_',x_name, '.csv');
    x=csvread(title_x,1,1);
    %x = importdata(title_x).data;
    title_y = strcat('02return_',y_name, '.csv');
    %y = importdata(title_y).data;
    y=csvread(title_y,1,1);
   [n,p] = size(x);   q = size(y,2);  
   [sA,sB] = gauss_data_ind(x,y);
   re1 = [n, p, q, sA, sB(1)];
   %a_1=num2cell(re1);
   %a2 = [x_name, y_name,a_1];
   %res_string(i,:)=a2;
   res(i,:)=re1 ;
   fprintf('the %d time \n', i)
end

title1 = ['02ind', '.csv'];
writetable(table(res),title1);


name = ["CS"  "CD" "CSt"  "Eng"  "Fin" "HC"  "Ind"  "IT"  "Mat"  "RE"  "Uti"];
com_my = nchoosek(name,2);
rep_num = size(com_my,1);
res = NaN(rep_num, 7);
%res_string = string(res);

for i=1:rep_num
    rng(12345);
    x_name = com_my(i,1);
    y_name = com_my(i,2);
    index1 = find(name == x_name);
    index2 = find(name == y_name);
    title_x = strcat('68return_',x_name, '.csv');
    x=csvread(title_x,1,1);
    %x = importdata(title_x).data;
    title_y = strcat('68return_',y_name, '.csv');
    %y = importdata(title_y).data;
    y=csvread(title_y,1,1);
   [n,p] = size(x);   q = size(y,2);  
   [sA,sB] = gauss_data_ind(x,y);
   re1 = [n, p, q, sA, sB(1)];
   %a_1=num2cell(re1);
   %a2 = [x_name, y_name,a_1];
   %%res_string(i,:)=a2;
   res(i,:)=re1 ;
    fprintf('the %d time \n', i)
end


title1 = ['68ind', '.csv'];
writetable(table(res),title1);
%xlswrite(title1,res_string) 


delete(gcp('nocreate'));