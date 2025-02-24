function [ x,y,z ] = data1_ex7( n,p,q,m,dep)
z = unifrnd (-1, 1, [n,m]);
x = randn(n, p);
y = randn(n,q);
rho = dep*0.1;
for i=1:m
    eps_y1 = unifrnd(-0.25, 0.25, [n, 48]);
    eps_y = sum(eps_y1, 2);
    y(:,i) = z(:,i) + 0.25 * z(:,i).^2 +  eps_y;
    beta = rho/(2*sqrt(1-rho^2));
    eps_x = randn(n,1);
    x(:,i)=5* beta .* y(:,i) + z(:,i) + eps_x;
end
for i=(p-m)
    x(:,(i+m)) = x(:,(i+m)) + 5*beta * sin( 2*y(:,(i+m)) );
end

end

