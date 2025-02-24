function [x,y] = data1_ex2(n,p,q,dep)
d = [floor(p/30),floor(p/25),floor(p/20),floor(p/15),floor(p/10)];
if(dep==0)
    x = 0.2* trnd(1, [n, p]);
    y = 0.2* trnd(1, [n, q]);
else
    x = 0.2* trnd(1, [n, p]);
    y = 0.2* trnd(1, [n, q]);
    d1 = d(dep);
    w = randn(n,1);
    x(:,1:d1)= x(:,1:d1) + w;
    y(:,1:d1)= y(:,1:d1) + w;
end
end

