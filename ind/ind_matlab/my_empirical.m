function y = my_empirical(x)
 n = length(x);
 [aa,bb]= ecdf(x);
 cc=zeros(n,2);
 cc(:,1)=aa(2:(n+1));
 cc(:,2)=bb(2:(n+1));
 y=zeros(n,1);
 for i=1:n
     ind = find(x==cc(i,2));
     y(ind)=cc(i,1);
 end
end

