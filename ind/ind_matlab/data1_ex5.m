function [x,y] = data1_ex5(n,p,q,dep)
d = [4,6,8,10,12];
if(dep==0)
    x = randn(n,p);
    y = randn(n,q);
else
    d1 = d(dep);
    Delta = zeros((p+q));
    comp = randsample(1:(p+q), d1);
    for i=1:(d1/2)
        index = (i-1)*2+1;
        a = unifrnd(0, 1, 1);
        Delta(comp(index),comp(index+1))= a;
        Delta(comp(index+1),comp(index))= a;
    end
    Id = eye(p+q);
    eival = min(eig(Id+Delta));
    delta = (-eival + 0.05) * (eival <=0);
    R_star = Id + Delta + delta* Id;
    z = mvnrnd(repelem(0,(p+q)),R_star,n);
    index1 = linspace(1,(d1-1),d1/2);
    index2 = linspace(2,d1,d1/2);
    x = zeros(n,p);
    y = zeros(n,q);
    x(:, 1:(d1/2)) = z(:, comp(index1));
    y(:, 1:(d1/2)) = z(:, comp(index2));
    %z1=z;
    z(:, comp)=[];
    x(:, (d1/2 +1): p) = z(:, 1:(p-d1/2));
    y(:, (d1/2 +1): q) = z(:, (p-d1/2+1):(p+q-d1));
end
for i=1:n
    for j=1:p
        x(i,j) = sign(x(i,j)) * abs(x(i,j))^(1/3);
        y(i,j) = sign(y(i,j)) * abs(y(i,j))^(1/3);
    end
end
end

