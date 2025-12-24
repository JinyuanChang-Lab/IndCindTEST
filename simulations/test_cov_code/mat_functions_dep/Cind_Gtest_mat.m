function out = Cind_Gtest_mat(x, y, z, alpha, option, N, seed)
% The proposed conditional independence test based on linear regressions with multiplier bootstrap
% Input:
%   x, y:   data matrices with the same number of rows (sample size n)
%   alpha:  significance level in (0, 1)
%   seed:   RNG seed (optional)
%   option: "Rademacher" | "Gaussian" | "Mammen" | "all"
%     "Rademacher" uses Rademacher multipliers;
%     "Gaussian" uses Gaussian multipliers;
%     "Mammen" uses Mammen's multipliers;
%     "all" computes results for all three types of multipliers
%   N: number of bootstrap replications, default 5000

% Output:
%   option:   the chosen option
%   reject:   1: reject the null hypothesis; 0: do not reject the null hypothesis
%   p_value:  p value
%   test_sta: test statistic
%   cv:       critical value

%%
% ensure x has the larger column dimension
if size(x, 2) < size(y, 2)
    tmp = x; x = y; y = tmp;
end

[n, p] = size(x);
q = size(y,2);
m = size(z,2);

% ---- coordinatewise Gaussianization ----
x_norm1 = data_transform_new(x);
y_norm1 = data_transform_new(y);
z_norm1 = data_transform_new(z);


%% lasso regression
lambda_min_u =NaN(1, p);
lambda_min_v =NaN(1, q);
a_hat =NaN(m,p);
b_hat =NaN(m,q);

% regress each column of x_norm1 on z_norm1 
parfor i= 1 : p
    rng(i*1000+123);
    fit1 = cvglmnet(z_norm1, x_norm1(:, i));       
    lambda_min_u(1,i) = fit1.lambda_min;        % select lambda that minimizes the CV error
    % calculate the coefficient with the selected lambda
    a_hat(:, i) = fit1.glmnet_fit.beta(:, find(fit1.lambda==fit1.lambda_min));
end

% regress each column of y_norm1 on z_norm1
parfor i = 1:q   
    rng(i*1000+123);
    fit2 = cvglmnet(z_norm1, y_norm1(:, i));        
    lambda_min_v(1,i) = fit2.lambda_min;         % select lambda that minimizes the CV error
    % calculate the coefficient with the selected lambda
    b_hat(:, i) = fit2.glmnet_fit.beta(:, find(fit2.lambda==fit2.lambda_min));
end

% calculate residuals
eps = x_norm1 - z_norm1 * a_hat;
delta = y_norm1 - z_norm1 * b_hat;


x_norm = eps;
y_norm = delta;

% ---- calculate the test statistic ----
test_sta = sqrt(n) * max(max(abs((x_norm.' * y_norm) / n)));

% --- calculate the critical values ---
% generate multipliers (N-by-n matrix)  
Mr = []; Mg = []; Mm = [];
needR = strcmpi(option,'Rademacher') || strcmpi(option,'all');
needG = strcmpi(option,'Gaussian')   || strcmpi(option,'all');
needM = strcmpi(option,'Mammen')     || strcmpi(option,'all');
if needR
    rng(seed);
    Mr = randsrc(N, n);                        % Rademacher multiplier 
end
if needG
    rng(seed);
    Mg = randn(N,n);                            % Gaussian multiplier 
end
if needM
    rng(seed);
    Mm = reshape(rmammen(N*n), N, n);            % Mammen's multiplier
end

% We need to split X into column-wise blocks to reduce memory usage when p and q are very large. tt is the block size (minimum 100)
tt = min(100, p);
B  = ceil(p / tt);           % the number of blocks
if needR, mxR = zeros(N, B); end
if needG, mxG = zeros(N, B); end
if needM, mxM = zeros(N, B); end

for b = 1:B
  i1 = (b-1)*tt + 1;
  i2 = min(b*tt, p);
  Xb = x_norm(:, i1:i2);                     
  tb = size(Xb,2);

  dat = reshape(repmat(Xb, q, 1), n, tb*q) .* repmat(y_norm, 1, tb);
  dat = dat - mean(dat,1);                  % Subtract the column mean from each column of dat

  if needR
    mxR(:,b) = max(abs(Mr * dat ./ sqrt(n)), [], 2);
  end
  if needG
    mxG(:,b) = max(abs(Mg * dat ./ sqrt(n)), [], 2);
  end
  if needM
    mxM(:,b) = max(abs(Mm * dat ./ sqrt(n)), [], 2);
  end
end
      
out = struct();


if strcmpi(option,'Rademacher') || strcmpi(option,'all')
  rowMax = max(mxR, [], 2);
  cvR    = sort(rowMax, 'descend'); 
  % critical value
  cvR1    = cvR(floor(N*alpha));
  % outputs: option, reject, tast statistics, p value
  outR   = pack('Rademacher', test_sta, cvR1, mean(rowMax > test_sta));
end
if strcmpi(option,'Gaussian') || strcmpi(option,'all')
  rowMax = max(mxG, [], 2);
  cvG    = sort(rowMax, 'descend'); 
  % critical value
  cvG1    = cvG(floor(N*alpha));
  % outputs: option, reject, tast statistics, p value
  outG   = pack('Gaussian', test_sta, cvG1, mean(rowMax > test_sta));
end
if strcmpi(option,'Mammen') || strcmpi(option,'all')
  rowMax = max(mxM, [], 2);
  cvM    = sort(rowMax, 'descend'); 
  % critical value
  cvM1    = cvM(floor(N*alpha));
  % outputs: option, reject, tast statistics, p value
  outM   = pack('Mammen', test_sta, cvM1, mean(rowMax > test_sta));
end

switch lower(option)
  case 'rademacher', out = outR;
  case 'gaussian',   out = outG;
  case 'mammen',     out = outM;
  otherwise
    out.Rademacher = outR;
    out.Gaussian   = outG;
    out.Mammen     = outM;
end
end



%% ---------- store the results ----------
function s = pack(kind, ts, cv, pval)
s = struct('option', kind, ...
           'reject', double(ts > cv), ...     % CV rule; use (pval<=alpha) if you prefer
           'p_value', pval, ...
           'test_sta', ts, ...
           'cv', cv);
end


%% ---------- generate mammen's multipliers ----------
function M_mm = rmammen(n)
  p = rand(n,1);                             % Generate n uniform random numbers in (0,1)
  reject = (sqrt(5)+1)/(2*sqrt(5));          % Threshold value for rejection  
  choice1 = -(sqrt(5)-1)/2;                   
  choice2 = (sqrt(5)+1)/2;                    
  M_mm = p;                                   
  for i = 1: n                                
      if (p(i) > reject)                     % If random number exceeds threshold
          M_mm(i) = choice2;                 % Assign choice2  
      else
          M_mm(i) = choice1;                 % Otherwise assign choice1  
      end
  end
end


%%  Define function that computes empirical CDF values for vector x
function y = my_empirical(x)  
 n = length(x);                       % Get the number of elements in x
 y=zeros(n,1);                         
 for i=1:n                             
     z=zeros(n,1);                     
     for j=1:n                         
     z(j) = x(j) <= x(i);             % Indicator: 1 if x(j) <= x(i), otherwise 0
     end                               
     y(i)=mean(z);                    % Empirical CDF at x(i): proportion of x(j) not greater than x(i)
 end                                   
end                                    

%%  coordinatewise Gaussianization
function [ x_t2 ] = data_transform_new(x)   
n = size(x,1);                             % the number of rows (observations)
p = size(x,2);                             % the number of columns (variables)
x_t2 = zeros(size(x));                     % Initialize output matrix with zeros of same size as x

for j=1 : p                                % Loop over each column
    x_t2(:,j) = my_empirical(x(:,j)) * n / (n+1);   % Calculate the empirical CDF, scaled by n/(n+1)
   
    for i=1 : n                                % Loop over each row in column j
        if (x_t2(i,j)~=0 && x_t2(i,j) ~=1 )    % If value is not 0 and not 1
           x_t2(i,j) = norminv(x_t2(i,j));     % standard normal quantile at x_t2(i, j)
        else
           x_t2(i,j) = 0;                      % Otherwise set value to 0
        end
    end
end

end                                         
