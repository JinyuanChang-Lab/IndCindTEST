function [ x_t2 ] = data_transform_new(x)
n = size(x,1);
p = size(x,2);
x_t2 = zeros(size(x));

 

for j=1 : p
    x_t2(:,j) = my_empirical(x(:,j)) *n/ (n+1) ;
   
    for i=1 : n
        if (x_t2(i,j)~=0 && x_t2(i,j) ~=1 )
           x_t2(i,j) = norminv(x_t2(i,j));
        else
           x_t2(i,j) = 0;
        end
    end
end

end



