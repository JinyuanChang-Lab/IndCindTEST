function [ x_t2 ] = data_transform( x,k )
n = size(x,1);
p = size(x,2);
M = sqrt(log(n)*k);
x_inf1 = max(normcdf(-M), 0);
x_sup1 = min(normcdf(M), 1);
x_t2 = zeros(size(x));

x_t2 = reshape(my_empirical(reshape(x, 1, n*p)), n, p) -1/n;
%%%%mimus 1/n
for j=1 : p
    %x_t2(:,j) = my_empirical(x(:,j)) - 1/n ;
    
    for i=1 : n
        if (x_t2(i,j) > x_inf1 &&  x_t2(i,j) <= x_sup1 &&  x_t2(i,j)~=0 && x_t2(i,j) ~=1 )
           x_t2(i,j) = norminv(x_t2(i,j));
        else
           x_t2(i,j) = 0;
        end
    end
end

end



