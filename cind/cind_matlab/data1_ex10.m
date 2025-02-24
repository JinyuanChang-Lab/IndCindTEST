function [ x,y,z ] = data1_ex10( n,p,q,m,dep)
z = randn(n,m);
w = randn(n,p);
v = randn(n,q);
dep=floor(dep);
x = randn(n,p); 
y=randn(n, q);
if (dep==0)
    for i= 1:(floor(p/4))
        zz = mean(z(:,1:m), 2);
        x(:,i) = tanh(zz + w(:,i));
        y(:,i) = (zz+v(:,i)).^3;
    end
else
    for i= 1:(floor(p/4))
        zz = mean(z(:,1:m), 2);
        x(:,i) = tanh(zz + w(:,i));
        y(:,i) = (zz + v(:,i)).^3;
    end
    for i= 1:dep
        eps = 3*randn(n, 1);
        zz = mean(z(:,1:m), 2);
        x(:,i) = tanh(zz + w(:,i)+ eps);
        y(:,i) = (zz + v(:,i) + eps ).^3;
    end

    
end

end
