function [x,y] = data1_ex3(n,p,q,dep)
d = [floor(p/20), floor(p/10)];
if(dep==0)
    x = unifrnd(0, 2*pi, [n, p]);
    y = unifrnd(0, 2*pi, [n, q]);
else
    x = unifrnd(0, 2*pi, [n, p]);
    y = unifrnd(0, 2*pi, [n, q]);
    d1 = d(dep);
    w = unifrnd(0, 2*pi, [n, d1]);
    x(:,1:d1)= (sin(w)).^2 ;
    y(:,1:d1)= (cos(w)).^2;
end
end

