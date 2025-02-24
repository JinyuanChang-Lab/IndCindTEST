function [x,y] = data1_ex1(n,p,q,dep)
d = [floor(p/20),floor(p/10)];
if(dep==0)
    x = trnd(1, [n, p]);
    y = trnd(1, [n, q]);
else
    x = trnd(1, [n, p]);
    y = trnd(1, [n, q]);
    d1 = d(dep);
    y(:,1:d1)= exp(x(:,1:d1));
end
end

