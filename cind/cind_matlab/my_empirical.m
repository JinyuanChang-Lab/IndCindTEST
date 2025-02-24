function y = my_empirical(x)
 n = length(x);
 y=zeros(n,1);
 for i=1:n
     z=zeros(n,1);
     for j=1:n
     z(j) = x(j) <= x(i); 
     end
     y(i)=mean(z);
 end
end

