%% ex1
clear
%%% Set the working directory to the current folder so that the glmnet
%%% package can be accessed properly
addpath(pwd); addpath(genpath(pwd));
rehash toolboxcache; savepath;
%choose cores
%parpool(8);


 


 
alpha=0.05;
%%%%%number
n_my = 2000;



 
% Scenario notes:
%   'indep_all'                 : X and Y independent; components within X and within Y are independent (S1 in Section R.3 of the supplementary material) 
%   'indep_vecs_corr_within'    : X and Y independent; components within X and/or Y are correlated (S2 in Section R.3 of the supplementary material) 
%   'corr_between_indep_within' : X and Y dependent; components within each are independent (S3 in Section R.3 of the supplementary material) 
%   'corr_all'                  : X and Y dependent; components within each also correlated (S4 in Section R.3 of the supplementary material) 
%   'dence'                     : X and Y dependent for each compnent (S5 in Section R.3 of the supplementary material) 
 
 

res_my = NaN(n_my,  1);
res_my_Schott = NaN(n_my,  1);
res_my_W = NaN(n_my,  1);
res_my_TBw = NaN(n_my,  1);
res_my_TBlh = NaN(n_my,  1);
res_my_TBbnp = NaN(n_my,  1);

 
 
p_values = [20, 40,  80, 200];  
n_values = [100, 200];  
m=5;


%% ---------- Settings S1-S4 ----------
scenarios = {'indep_all', 'indep_vecs_corr_within', 'corr_between_indep_within', 'corr_all'};  
for n_idx = 1:length(n_values)
    n = n_values(n_idx);
    fprintf('-----------  n %d ----------- \n', n)
    for scenario_idx = 1:length(scenarios)
        scenario = scenarios{scenario_idx};   
         
        for idx = 1:length(p_values)
            p1 = p_values(idx);   
            p2 = p1;              
            fprintf('-----------  n %d ----------- \n', n);
            fprintf('-----------  p1 %d ----------- \n', p1);
            fprintf('-----------  scenario %s ----------- \n', scenario);
            parfor i=1 : n_my
                %tic
                rng(i*1000+123*scenario_idx + idx);
                %generate data 
                [x,y,z,a1,a2,SigE] = make_XY_once(n,p1,p1,m,scenario);
                % calcutale the residuals of x
                A1    = z \ x;       
                R1    = x - z * A1;
                % calcutale the residuals of y
                A2    = z \ y;       
                R2    = y - z * A2;
                
                % the proposed method
                res1 = Cind_Gtest_mat(x, y, z, alpha, 'Rademacher', 5000, n*1000+123);
                res_my(i,:) = res1.reject;


                %residual of X and Y
                X = [R1, R2];

                try
                    % (Schott) in Bao et al. (2017) 
                    Schott = Schott_statistics(X, p1, p2);
                    if ~isreal(Schott) || isnan(Schott)
                        warning('Schott result is complex or NaN at iteration %d', i);
                        res_my_Schott(i, 1) = NaN;  % indicates failure
                    else
                        res_my_Schott(i, 1) = (Schott > 1.6448);
                    end

                    % (LRT) in Jiang and Yang (2013)
                    W = Wilk_statistic(X, p1, p2);
                    if ~isreal(W) || isnan(W)
                        warning('Wilks result is complex or NaN at iteration %d', i);
                        res_my_W(i, 1) = NaN;  % indicates failure 
                    else
                        res_my_W(i, 1) = (W > 1.6448);
                    end

                    %  TW, TLH,TBNP in Bodnar et al. (2019) 
                    [W_w, W_lh, W_bnp] = TB_statistics(X, p1, p2);

                    % Check for complex or NaN values in the TB statistics
                    if ~isreal(W_w) || isnan(W_w)
                        warning('Wilks result is complex or NaN at iteration %d', i);
                        res_my_TBw(i, 1) = NaN;  % indicates failure 
                    else
                        res_my_TBw(i, 1) = (W_w > 1.6448);
                    end
                    if ~isreal(W_lh) || isnan(W_lh)
                        warning('Wilks result is complex or NaN at iteration %d', i);
                        res_my_TBlh(i, 1) = NaN;  % indicates failure 
                    else
                        res_my_TBlh(i, 1) = (W_lh > 1.6448);
                    end

                    if ~isreal(W_bnp) || isnan(W_bnp)
                        warning('Wilks result is complex or NaN at iteration %d', i);
                        res_my_TBbnp(i, 1) = NaN;  % indicates failure 
                    else
                        res_my_TBbnp(i, 1) = (W_bnp > 1.6448);
                    end
                catch exception
                    % error
                    warning('Error occurred at iteration %d: %s', i, exception.message);
                end
             %toc
            %writetable(table(res_my),title);    
                %fprintf('the %d time \n', i)
            end


            %% Compute the proportion of rejections of the null across n_my replications
            result = mean([res_my(:,1), res_my_Schott,res_my_W, res_my_TBw,res_my_TBlh,res_my_TBbnp], 1);
            %title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
            title1 = sprintf('%s_n_%d_p_%d_q_%d%d.csv', scenario, n, p1, p2);
            title2 = [title1, '.csv'];
            writetable(table(result),title2);

       
        end
    end
end 
 


%% ---------- Setting S5 ----------
for n_idx = 1:length(n_values)
    n = n_values(n_idx);
    fprintf('-----------  n %d ----------- \n', n)
    scenario = 'dence';
    for idx = 1:length(p_values)
        p1 = p_values(idx);   
        p2 = p1;              
        fprintf('-----------  n %d ----------- \n', n);
        fprintf('-----------  p1 %d ----------- \n', p1);
        fprintf('-----------  scenario %s ----------- \n', scenario);
        parfor i=1 : n_my
            %tic
            rng(i*1000+123);

            [x,y,z,a1,a2] = make_XY_once_dence(n,p1,p2,m);
            A1    = z \ x;       
            R1    = x - z * A1;
            A2    = z \ y;       
            R2    = y - z * A2;


            % the proposed method
            res1 = Cind_Gtest_mat(x, y, z, alpha, 'Rademacher', 5000, n*1000+123);
            res_my(i,:) = res1.reject;  


            %residual of X and Y
            X = [R1, R2];

            try
                % (Schott) in Bao et al. (2017) 
                Schott = Schott_statistics(X, p1, p2);
                if ~isreal(Schott) || isnan(Schott)
                    warning('Schott result is complex or NaN at iteration %d', i);
                    res_my_Schott(i, 1) = NaN;    % indicates failure
                else
                    res_my_Schott(i, 1) = (Schott > 1.6448);
                end

                % (LRT) in Jiang and Yang (2013)
                W = Wilk_statistic(X, p1, p2);
                if ~isreal(W) || isnan(W)
                    warning('Wilks result is complex or NaN at iteration %d', i);
                    res_my_W(i, 1) = NaN;   % indicates failure
                else
                    res_my_W(i, 1) = (W > 1.6448);
                end

                %  TW, TLH,TBNP in Bodnar et al. (2019) 
                [W_w, W_lh, W_bnp] = TB_statistics(X, p1, p2);

                 
                if ~isreal(W_w) || isnan(W_w)
                    warning('Wilks result is complex or NaN at iteration %d', i);
                    res_my_TBw(i, 1) = NaN;   % indicates failure
                else
                    res_my_TBw(i, 1) = (W_w > 1.6448);
                end
                if ~isreal(W_lh) || isnan(W_lh)
                    warning('Wilks result is complex or NaN at iteration %d', i);
                    res_my_TBlh(i, 1) = NaN;   % indicates failure
                else
                    res_my_TBlh(i, 1) = (W_lh > 1.6448);
                end

                if ~isreal(W_bnp) || isnan(W_bnp)
                    warning('Wilks result is complex or NaN at iteration %d', i);
                    res_my_TBbnp(i, 1) = NaN;   % indicates failure
                else
                    res_my_TBbnp(i, 1) = (W_bnp > 1.6448);
                end
            catch exception
                % error
                warning('Error occurred at iteration %d: %s', i, exception.message);
            end
         %toc
        %writetable(table(res_my),title);    
            %fprintf('the %d time \n', i)
        end

        result = mean([res_my(:,1), res_my_Schott,res_my_W, res_my_TBw,res_my_TBlh,res_my_TBbnp], 1);
        %title1 = ['ex',num2str(ex) ,'_n_', num2str(n), '_p_', num2str(p), '_dep_', num2str(dep), '.csv'];
        title1 = sprintf('%s_n_%d_p_%d_q_%d%d.csv', scenario, n, p1, p2);
        title2 = [title1, '.csv'];
        writetable(table(result),title2);


    end
end 



delete(gcp('nocreate'));