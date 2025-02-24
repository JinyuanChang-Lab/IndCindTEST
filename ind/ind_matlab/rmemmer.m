function M_mm = rmemmer(n)
  p = rand(n,1);
  reject = (sqrt(5)+1)/(2*sqrt(5));
  choice1 = -(sqrt(5)-1)/2;
  choice2 = (sqrt(5)+1)/2;
  M_mm = p;
  for i = 1: n
      if (p(i) > reject)
          M_mm(i) = choice2;
      else
          M_mm(i) = choice1;
      end
  end
end

