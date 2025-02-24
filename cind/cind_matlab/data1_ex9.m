function [ x,y,z ] = data1_ex9( n,p,q,m,dep)
z = randn(n,m);
w = randn(n,p);
v = randn(n,q);
x = randn(n,p);
y = randn(n, q);
for i =1:m
    a = 0.5 .*(z(:, i).^3./7 + z(:, i)./2);
    b = (z(:,i).^3./2 + z(:,i))./3;
    a1 = a+ tanh(w(:,i));
    x(:,i) = a1 + a1.^3./3;
    b1 = b + v(:,i);
    y(:,i) = b1 + tanh(b1/3); 
end
if dep~=0
    for i= 1:dep
        eps = randn(n, 1);
        x(:,i)= 0.5*x(:,i) + 3*cosh(eps);
        y(:,i)= 0.5*y(:,i) + 3*cosh(eps.^2);
    end
end

end


