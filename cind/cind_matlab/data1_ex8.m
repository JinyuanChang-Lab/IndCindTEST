function [ x,y,z ] = data1_ex8( n,p,q,m,dep)
z = randn(n,m);
w = randn(n,p);
v = randn(n,q);
x = randn(n,p); 
y = randn(n, q);
for i =1:m
    a = 0.7 .*(z(:, i).^3./5 + z(:, i)./2) + tanh((w(:,i)));
    b = (z(:,i).^3./4 + z(:,i))./3 + v(:,i);
    x(:,i) = a + a.^3./3 + tanh(a./3)./2;
    y(:,i) = (b + tanh(b./3)).^3; 
end

if dep ~=0
    for i= 1:dep
        eps = 3*trnd(1, [n, 1]);
        x(:,i)= x(:,i) + eps;
        y(:,i)= y(:,i) + eps;
    end
end

end


