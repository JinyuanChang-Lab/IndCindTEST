"""
The function `Cind_Gtest_py` implements the proposed conditional independence test based on nonparametric regressions.
Args:
    setting: (1) "S-selectn3":  splitting sample with selecting n_3^{opt}  (see Algorithm, Section 4.1).
             (2) "S-givenn3":   splitting sample but without selecting n_3^{opt} (see Remark 4, Section 4.1).
             (3) "full-sample": without sample-splitting (see Remark 4, Section 4.1).
    option_l: Rademacher / Gaussian / Mammen
""" 

import numpy as np
import torch
import torch.nn as nn

from torch.utils.data import TensorDataset, DataLoader, random_split
import math

from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF

import logging
import matplotlib.pyplot as plt


# the proposed functions
import non_connect_nn


def init_logging(log_file):
    """
    Define the log function
    """

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
     
    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler.setFormatter(formatter)
    
     
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    
    if not logger.handlers:
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
    
    return logger



def generate_one_minus_one_distribution(a, b):
    """
    generate Radamacher multiplier
    """
    return 2 * torch.randint(0, 2, (a, b)) - 1
    


def generate_mammen(a, b):
    """
    generate Mammen's multiplier
    """
    pv = torch.rand(a, b)
    reject = (np.sqrt(5) + 1) / (2 * np.sqrt(5))
    choice1 = -(np.sqrt(5) - 1) / 2
    choice2 = (np.sqrt(5) + 1) / 2

    mask = pv > reject
    M_mm = np.where(mask, choice2, choice1)

    return torch.tensor(M_mm)



def data_transform(x_1,x_2):
    """
    coordinatewise Gaussianization
    """
    # x_1 is an n_1-by-p matrix, x_2 is an n_2-by-p matrix
    n_1 = x_1.shape[0]                                   # the number of rows of x_1
    dim_x = x_2.shape                                    # the dimension of x_2
    pp = dim_x[1]                                        # the number of columns of x_2
    y_1 = torch.zeros(dim_x)                             # generate an n_2-by-p matrix

    # calculate the empirical CDF of each column of x_1
    ecdf = [ECDF(x_1[:, i]) for i in range(x_1.shape[1])]   

    # calculate the truncated empirical CDF given in (11) in Section 4.1.1 of the main paper
    for j in range(0,pp):
        y_1[:, j] = torch.tensor(ecdf[j](x_2[:, j]), dtype=torch.float32)  # evaluate the empirical CDF of x_1[:,j] at x_2[:,j]
        y_1[:, j] = torch.clamp(y_1[:, j], 1/n_1, (1-1/n_1))               # truncate the ECDF into (1/n1, 1-1/n1) to avoid taking values 0 and 1
        # y_1[:, j]= my_empirical(x_1[:, j], x_2[:, j]).T                   
        y_1[:, j] = torch.tensor(norm.ppf(y_1[:,j], 0, 1))                 # standard normal quantile at y_1[:,j]

    return y_1         
    


def threshold(mat, thresh):
    """
    Restrict all elements of mat to the interval [-thresh, thresh]
    """
    mat = torch.clamp(mat, min=-thresh, max=thresh)
    return mat



def trainer(setting, exp_v, res_v, n, batchsize, lr, n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, device, thresh=None, drop_last1=True, plot=False):
    """
    Train a two-hidden-layer network with a fixed row-wise split and early stopping.

    Procedure:
        1) Use early stopping and ReduceLROnPlateau on validation loss.
        2) Return the fitted model and residuals. 

    Args: 
        exp_v: Covariate matrix of shape (n, m).
        res_v: Response vector of shape (n, 1).
    """
     
    # ----------------------------
    # 1) Setup: model, optimizer, scheduler
    # ----------------------------
    criterion = nn.MSELoss()
    model = non_connect_nn.AdaptiveNet(input_features, hidden_features1, hidden_features2, output_features).to(device)
    
    # optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    optimizer = torch.optim.Adam(model.parameters(),
                                   lr=lr,
                                   weight_decay=0.)
    
    # If validation loss does not improve for 5 epochs, reduce the learning rate by half (not below 1e-4).
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=5, factor=0.5, min_lr=1e-4)
     

    # ----------------------------
    # 2) Data split 
    # ----------------------------
    if setting in ("S-selectn3", "S-givenn3"):
        n_1 = int(n/3)
        n_2 = int(n/2)
        # train on [n1:n1+n2]; validate on the remaining data ([0:n1] and [n1+n2:n]).
        train_dataset = TensorDataset(exp_v[n_1:(n_1+n_2),:], res_v[n_1:(n_1+n_2)])    # training set
    
        n1_ind = torch.arange(0, n_1)
        n2_ind = torch.arange(n_1+n_2, n)
        all_indices = torch.cat([n1_ind, n2_ind])
        val_dataset = TensorDataset(exp_v[all_indices,:], res_v[all_indices])     # validation set
    else: # for full-sample
        # train on [0:0.8*n]; validate on the remaining data ([0.8*n:n]).
        train_dataset = TensorDataset(exp_v[0:int(n*0.8),:], res_v[0:int(n*0.8)])
        all_indices = torch.arange(int(n*0.8), n)
        val_dataset = TensorDataset(exp_v[all_indices,:], res_v[all_indices])
     
    train_loader = DataLoader(train_dataset, batch_size=batchsize, shuffle=True,drop_last=drop_last1)
    val_loader = DataLoader(val_dataset, batch_size=batchsize, shuffle=False)


    # ----------------------------
    # 3) Training loop with early stopping on validation loss
    # ----------------------------
    run_loss = 0   
    min_val_loss = float('inf')   
    no_improvement_epochs = 0     # no-improvement counter

    val_losses = []   
    run_losses = []   
    for epoch in range(n_epochs):
        model.train()   
        run_loss1 = run_loss   
        run_loss = 0   
        # Train the network using data from train_loader 
        for inputs, labels in train_loader:
            inputs = inputs.to(device)   
            labels = labels.to(device)   
            optimizer.zero_grad()   
            outputs = model(inputs)   
            loss = criterion(outputs.view(-1), labels)   
            loss.backward()   
            optimizer.step()   
            run_loss += loss.item()   

        run_loss = run_loss / len(train_loader)   
        run_losses.append(run_loss)   
         

        # ----- validation -----
        model.eval()   
        val_loss = 0   
        with torch.no_grad():  # no gradients during validation
            for val_inputs, val_labels in val_loader:
                val_inputs = val_inputs.to(device)   
                val_labels = val_labels.to(device)   
                val_outputs = model(val_inputs)   
                val_loss += criterion(val_outputs.view(-1), val_labels).item()   
        val_loss = val_loss / len(val_loader)   
        val_losses.append(val_loss)   
        scheduler.step(val_loss)   

        if val_loss < min_val_loss:   
            min_val_loss = val_loss                                                      # If the current validation loss improves, update min_val_loss; Otherwise, increase the no-improvement counter
            no_improvement_epochs = 0   
            best_model_state = {k: v.clone() for k, v in model.state_dict().items()}     # save best weights
            # torch.save(model_all.state_dict(), 'best_model.pth')   
        else:
            no_improvement_epochs += 1                             

        # Early stopping rules: stop if no improvement for 'patience' epochs or if min_val_loss is very small
        if no_improvement_epochs >= patience or min_val_loss < 1e-3:
            break   

    # ----------------------------
    # 4) Calculate the residuals   
    # ----------------------------
    model.load_state_dict(best_model_state)   
    if thresh:   
        with torch.no_grad():
            model.eval()
            fit = threshold(model.forward(exp_v.to(device)), thresh)   
    else:
        with torch.no_grad():
            model.eval()
            fit = model(exp_v.to(device))   
    
    # Calculate the residuals  
    error = res_v.to(device) - fit.squeeze()   

    # ----------------------------
    # 5) Plot learning curves
    # ----------------------------
    if plot:   
        fig = plt.figure(figsize=(10, 4))
        plt.style.use('fivethirtyeight')
        plt.plot(run_losses, label='Training Loss', color='C0', linewidth=2)
        plt.plot(val_losses, label='Validation Loss', color='C1', linewidth=2)
        plt.yscale('log')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return {
        'model': model,         
        'error': error,         
        'losses': run_losses,   
        'val_losses': val_losses,   
    }


def gauss_fun_split(setting, option_l, x, y, z, n, alpha,batchsize, lr,n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features,  drop_last1, device):
    """
    implement the proposed testing procedure with sample-splitting (see Algorithm 1 and Remark 4 in Section 4.1) 
    """
    n_1 = int(n/3)
    n_2 = int(n/2)
    n_3 = n - n_1 -n_2
    n_tilde = torch.arange(int(n_3/3), n_3, 2)
    p = x.size(1)
    q = y.size(1)
    m = z.size(1)
    
    # coordinatewise Gaussianization
    x = data_transform(x[0:n_1,:], x)
    y = data_transform(y[0:n_1,:], y)
    z = data_transform(z[0:n_1,:], z)

    # Concatenate response
    rep_xy = torch.cat([x, y], dim=1)
    
    m1=max(p,q,m)   

    thresh = (math.log(m1*n))**(1/2) * math.log(n)

    # For every column in [x, y], fit a neural network using z as predictors
    results = [trainer(setting, z, rep_xy[:, i], n, batchsize,  lr, n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, device,thresh, drop_last1) for i in range(p + q)]
    
    # calculate residuals  
    error_all = torch.stack([res['error'] for res in results],dim=1)   # stack residuals from all (p+q) models into an n-by(p+q) matrix 
    eps_all = error_all[:, 0:p]           # get residuals of x on z
    delta_all = error_all[:, p:(p+q)]     # get residuals of y on z
    
    ell = n_tilde.shape[0]                # length of n_tilde
    N=5000                                # number of multiplier bootstrap replications
    N_1 = 500                             # number of data-driven replications to choose n_3^{opt}
    res1_red = torch.zeros(N_1, ell)      # generate an N_1-by-ell matrix
    res1_norm = torch.zeros(N_1, ell)     # generate an N_1-by-ell matrix
    res1_mm = torch.zeros(N_1, ell)       # generate an N_1-by-ell matrix

    #--------------------------------------------- setting  "S-givenn3" -------------------
    if setting == "S-selectn3":
        # ------------------------------------------------------------------------------------------------ 
        # Data-driven selection of n_3^{opt} via Algorithm 1 with Rademacher/Gaussian/Mammen's multipliers 
        # ------------------------------------------------------------------------------------------------
        for t in range(N_1):   
            f1 = (torch.randn(n, n) / (n**(1/2))).to(device)       # generate an n-by-n matrix with entries drawn from N(0, 1), scaled by n^{-1/2}
            f2 = (torch.randn(n, n) / (n**(1/2))).to(device)       # generate an n-by-n matrix with entries drawn from N(0, 1), scaled by n^{-1/2}
            f3 = (torch.randn(n, n) / (n**(1/2))).to(device)       # generate an n-by-n matrix with entries drawn from N(0, 1), scaled by n^{-1/2}
            w_hs_norm = f1 @ z.to(device)  # matrix product: f1 and z
            eps_norm = f2 @ eps_all        # matrix product: f2 and eps_all
            delta_norm = f3 @ delta_all    # matrix product: f3 and delta_all
            # generate date used in Alogrithm 1
            with torch.no_grad():    # no gradients
                f_h_hs = torch.stack([res['model'].forward(w_hs_norm.to(device)).squeeze() for res in results],dim=1)   # compute network outputs given input w_hs_norm 
            thresh1 = norm.ppf((1-1/n_1), 0, 1)
            ###generate hat{w}, hat{u} and hat{v} in the Algorithm 1
            w_hs_norm  = threshold(w_hs_norm, thresh1)
            u_hs_norm = f_h_hs[:, 0: p] + eps_norm
            u_hs_norm = threshold(u_hs_norm, thresh1)
            v_hs_norm = f_h_hs[:, p: (p+q)] + delta_norm
            v_hs_norm = threshold(v_hs_norm, thresh1)
            ## train neural networks with the data hat{w}, hat{u} and hat{v}
            if t<=0:
                rep_uv = torch.cat([u_hs_norm, v_hs_norm], dim=1)
                results2 = [trainer(setting, w_hs_norm, rep_uv[:, k], n, batchsize,  lr, n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, device,thresh, drop_last1) for k in range(p + q)]
             
            with torch.no_grad():    # no gradients
                f_hhm = torch.stack([res['model'].forward(w_hs_norm.to(device)).squeeze() for res in results2],dim=1)
            
            ## calculate the residuals
            error_all_u = torch.cat([u_hs_norm, v_hs_norm], dim=1) - f_hhm
    
            ### select rows n1+n2 : n from error_all_u  
            eps_hhs = error_all_u[(n_1+n_2):n, 0:p]
            delta_hhs = error_all_u[(n_1+n_2):n, p:(p+q)] 

            ### generate Radamacher multipliers (N-by-n_3 matrix)   
            M = generate_one_minus_one_distribution(N, n_3).float().to(device)
            ### generate Gaussian multipliers (N-by-n_3 matrix)
            M_norm = torch.randn(N, n_3).to(device) 
            ### generate Mamman's multipliers (N-by-n_3 matrix)
            M_mm = generate_mammen(N,n_3).float().to(device)
    
            # For each candidate n_tilde[i], compute bootstrap critical values and indicator of rejection 
            for i in range(ell):
                
                n_tilde1 = n_tilde[i]
                sqrt_n_tilde1 = torch.sqrt(torch.as_tensor(n_tilde1)).to(device)
                dat_new = torch.empty(N, p).to(device)      # generate an N-by-p matrix
                dat_new_n = torch.empty(N, p).to(device)    # generate an N-by-p matrix
                dat_new_mm = torch.empty(N, p).to(device)   # generate an N-by-p matrix
                # process eps_hhs columnwise to reduce memory usage
                for j in range(p):
                    Ai = eps_hhs[0:n_tilde1, j:j+1]
                    dat_i = Ai * delta_hhs[0:n_tilde1, :]   
                    dat_i = (dat_i - torch.mean(dat_i, dim=0, keepdim=True)).to(device)
                    if option_l in ("Rademacher", "all"):
                        ### Radamacher multiplier   
                        result_i = torch.max(abs(M[:, 0:n_tilde1] @ dat_i / sqrt_n_tilde1), axis=1)[0]   
                        dat_new[:, j] = result_i   
                    if option_l in ("Gaussian", "all"):    
                        ### Gaussian multiplier
                        result_i_n = torch.max(abs(M_norm[:, 0:n_tilde1] @ dat_i / sqrt_n_tilde1), axis=1)[0]   
                        dat_new_n[:, j] = result_i_n   
                    if option_l in ("Mammen", "all"):    
                        ### Mammen's multiplier
                        result_i_mm = torch.max(abs(M_mm[:, 0:n_tilde1] @ dat_i / sqrt_n_tilde1), axis=1)[0]   
                        dat_new_mm[:, j] = result_i_mm   
    
                # calculate the test statistic 
                test1 = np.sqrt(n_tilde1) * torch.max(abs(eps_hhs[0:n_tilde1, :].T @ delta_hhs[0:n_tilde1, :] / n_tilde1))
                    
                if option_l in ("Rademacher", "all"):        
                    dat_new = torch.max(dat_new, axis = 1).values         # calculate the row-wise maximum
                    # calculate the critical value   
                    cv_red = torch.kthvalue(dat_new, (N-int(N*alpha)+1)).values
                    # calculate the p value 
                    pval_red = (dat_new > test1).float().mean()
                    # Reject the null hypothesis (=1) if the test statistic > critical value
                    res1_red[t, i]  = (test1 >cv_red)
    
                if option_l in ("Gaussian", "all"):
                    dat_new_n = torch.max(dat_new_n, axis = 1).values     # calculate the row-wise maximum
                    # calculate the critical value  
                    cv_norm = torch.kthvalue(dat_new_n, (N-int(N*alpha)+1)).values
                    # Reject the null hypothesis (=1) if the test statistic > critical value
                    res1_norm[t, i] = (test1 >cv_norm)
                    
                if option_l in ("Mammen", "all"):   
                    dat_new_mm = torch.max(dat_new_mm, axis = 1).values   # calculate the row-wise maximum
                    # calculate the critical value with Mammen's multipliers
                    cv_mm = torch.kthvalue(dat_new_mm, (N-int(N*alpha)+1)).values
                    # Reject the null hypothesis (=1) if the test statistic > critical value
                    res1_mm[t, i]   = (test1 >cv_mm)
                
        # --------------------------------------------------------
        # Choose n_3^{opt} and calculate test statistics with it  
        # --------------------------------------------------------
        n_s = n_1 + n_2
        ####Rademacher  
        if option_l in ("Rademacher", "all"): 
            size_red  = torch.mean(res1_red, dim=0)    # empirical size = no. of rejects / N 
            index_red = int(max(torch.where(abs(size_red-alpha) == min(abs(size_red-alpha)))[0]))   
            n_tilde_red = int(n_tilde[index_red])      # choose n_3^{opt} (Rademacher) with empirical size closest to alpha
        
            eps_t1 = eps_all[n_s:(n_s+n_tilde_red), :]
            delta_t1 = delta_all[n_s:(n_s+n_tilde_red), :]  
            test_1 = torch.max(abs(eps_t1.T @ delta_t1 / np.sqrt(n_tilde_red) ))   # calculate test statistics with selected n_3^{opt} (Rademacher)
            ### generate Radamacher multipliers (N-by-n_tilde_red matrix) 
            M = generate_one_minus_one_distribution(N, n_tilde_red).float().to(device)
            
        ####Gaussian
        if option_l in ("Gaussian", "all"):
            size_norm  = torch.mean(res1_norm, dim=0)   # empirical size = no. of rejects / N
            index_norm = int(max(torch.where(abs(size_norm-alpha) == min(abs(size_norm-alpha)))[0]))
            n_tilde_norm = int(n_tilde[index_norm])     # choose n_3^{opt} (Gaussian) with empirical size closest to alpha
    
            eps_t2 = eps_all[n_s:(n_s+n_tilde_norm), :]
            delta_t2 = delta_all[n_s:(n_s+n_tilde_norm), :]
            test_2 =  torch.max(abs(eps_t2.T @ delta_t2 / np.sqrt(n_tilde_norm) ))  # calculate test statistics with selected n_3^{opt} (Gaussian)
            ### generate Gaussian multipliers (N-by-n_tilde_norm matrix)
            M_norm = torch.randn(N, n_tilde_norm).to(device)
         
        ####Mammen's    
        if option_l in ("Mammen", "all"):
            size_mm  = torch.mean(res1_mm, dim=0)   # empirical size = no. of rejects / N
            index_mm = int(max(torch.where(abs(size_mm-alpha) == min(abs(size_mm-alpha)))[0]))
            n_tilde_mm = int(n_tilde[index_mm])      # choose n_3^{opt} (Mammen's) with empirical size closest to alpha
        
            eps_t3 = eps_all[n_s:(n_s+n_tilde_mm), :]
            delta_t3 = delta_all[n_s:(n_s+n_tilde_mm), :]
            test_3 = torch.max(abs(eps_t3.T @ delta_t3 / np.sqrt(n_tilde_mm) )) # calculate test statistics with selected n_3^{opt} (Mammen's)
            ### generate Mammen's multipliers (N-by-n_tilde_mm matrix)
            M_mm = generate_mammen(N,n_tilde_mm).float().to(device)
        
    
        # ------------------------------
        # Calculate the critical values  
        # ------------------------------ 
        if option_l in ("Rademacher", "all"):
            n_tilde_red = torch.sqrt(torch.as_tensor(n_tilde_red)).to(device)   
        if option_l in ("Gaussian", "all"):
            n_tilde_norm = torch.sqrt(torch.as_tensor(n_tilde_norm)).to(device)  
        if option_l in ("Mammen", "all"):
            n_tilde_mm = torch.sqrt(torch.as_tensor(n_tilde_mm)).to(device)   
    
        dat_new1 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
        dat_new2 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
        dat_new3 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
    
        for i in range(p):
            ### Radamacher multiplier 
            if option_l in ("Rademacher", "all"): 
                dat_1 = eps_t1[:, i:i+1] * delta_t1   
                dat_1 = (dat_1 - torch.mean(dat_1, dim=0, keepdim=True)).to(device)
                dat_new1[:, i] = torch.max(abs(M @ dat_1 / n_tilde_red), axis=1)[0]  
            ### Gaussian multiplier
            if option_l in ("Gaussian", "all"):
                dat_2 = eps_t2[:, i:i+1] * delta_t2   
                dat_2 = (dat_2 - torch.mean(dat_2, dim=0, keepdim=True)).to(device)   
                dat_new2[:, i] = torch.max(abs(M_norm @ dat_2 / n_tilde_norm), axis=1)[0]   
            ### Mammen's multiplier
            if option_l in ("Mammen", "all"):
                dat_3 = eps_t3[:, i:i+1] * delta_t3   
                dat_3 = (dat_3 - torch.mean(dat_3, dim=0, keepdim=True)).to(device)   
                dat_new3[:, i] = torch.max(abs(M_mm @ dat_3 / n_tilde_mm), axis=1)[0]   
     
        if option_l in ("Rademacher", "all"):
            dat_new1 = torch.max(dat_new1, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Radamacher multiplier
            cv_red1 = torch.kthvalue(dat_new1, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_red1 = (dat_new1 > test_1).float().mean()
            rej_red1 = (test_1 > cv_red1)
        if option_l in ("Gaussian", "all"):
            dat_new2 = torch.max(dat_new2, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Gaussian multiplier
            cv_norm1 = torch.kthvalue(dat_new2, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_norm1 = (dat_new2 > test_2).float().mean()
            rej_norm1 = (test_2 > cv_norm1)
        if option_l in ("Mammen", "all"):    
            dat_new3 = torch.max(dat_new3, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Mammen's multiplier
            cv_mm1 = torch.kthvalue(dat_new3, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_mm1 = (dat_new3 > test_3).float().mean()
            rej_mm1 = (test_3 > cv_mm1)

    #--------------------------------------------- setting  "S-givenn3" -------------------
    if setting == "S-givenn3":
        n_s = n_1 + n_2
        ####Rademacher  
        if option_l in ("Rademacher", "all"): 
            n_tilde_red = n-n_s
            eps_t1 = eps_all[n_s:(n_s+n_tilde_red), :]
            delta_t1 = delta_all[n_s:(n_s+n_tilde_red), :]  
            test_1 = torch.max(abs(eps_t1.T @ delta_t1 / np.sqrt(n_tilde_red) ))   # calculate test statistics with selected n_3^{opt} (Rademacher)
            ### generate Radamacher multipliers (N-by-n_tilde_red matrix) 
            M = generate_one_minus_one_distribution(N, n_tilde_red).float().to(device)
            
        ####Gaussian
        if option_l in ("Gaussian", "all"):
            n_tilde_norm = n-n_s
            eps_t2 = eps_all[n_s:(n_s+n_tilde_norm), :]
            delta_t2 = delta_all[n_s:(n_s+n_tilde_norm), :]
            test_2 =  torch.max(abs(eps_t2.T @ delta_t2 / np.sqrt(n_tilde_norm) ))  # calculate test statistics with selected n_3^{opt} (Gaussian)
            ### generate Gaussian multipliers (N-by-n_tilde_norm matrix)
            M_norm = torch.randn(N, n_tilde_norm).to(device)
         
        ####Mammen's    
        if option_l in ("Mammen", "all"):
            n_tilde_mm = n-n_s
            eps_t3 = eps_all[n_s:(n_s+n_tilde_mm), :]
            delta_t3 = delta_all[n_s:(n_s+n_tilde_mm), :]
            test_3 = torch.max(abs(eps_t3.T @ delta_t3 / np.sqrt(n_tilde_mm) )) # calculate test statistics with selected n_3^{opt} (Mammen's)
            ### generate Mammen's multipliers (N-by-n_tilde_mm matrix)
            M_mm = generate_mammen(N,n_tilde_mm).float().to(device)
        
    
        # ------------------------------
        # Calculate the critical values  
        # ------------------------------ 
        if option_l in ("Rademacher", "all"):
            n_tilde_red = torch.sqrt(torch.as_tensor(n_tilde_red)).to(device)   
        if option_l in ("Gaussian", "all"):
            n_tilde_norm = torch.sqrt(torch.as_tensor(n_tilde_norm)).to(device)  
        if option_l in ("Mammen", "all"):
            n_tilde_mm = torch.sqrt(torch.as_tensor(n_tilde_mm)).to(device)   
    
        dat_new1 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
        dat_new2 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
        dat_new3 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
    
        for i in range(p):
            ### Radamacher multiplier 
            if option_l in ("Rademacher", "all"): 
                dat_1 = eps_t1[:, i:i+1] * delta_t1   
                dat_1 = (dat_1 - torch.mean(dat_1, dim=0, keepdim=True)).to(device)
                dat_new1[:, i] = torch.max(abs(M @ dat_1 / n_tilde_red), axis=1)[0]  
            ### Gaussian multiplier
            if option_l in ("Gaussian", "all"):
                dat_2 = eps_t2[:, i:i+1] * delta_t2   
                dat_2 = (dat_2 - torch.mean(dat_2, dim=0, keepdim=True)).to(device)   
                dat_new2[:, i] = torch.max(abs(M_norm @ dat_2 / n_tilde_norm), axis=1)[0]   
            ### Mammen's multiplier
            if option_l in ("Mammen", "all"):
                dat_3 = eps_t3[:, i:i+1] * delta_t3   
                dat_3 = (dat_3 - torch.mean(dat_3, dim=0, keepdim=True)).to(device)   
                dat_new3[:, i] = torch.max(abs(M_mm @ dat_3 / n_tilde_mm), axis=1)[0]   
     
        if option_l in ("Rademacher", "all"):
            dat_new1 = torch.max(dat_new1, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Radamacher multiplier
            cv_red1 = torch.kthvalue(dat_new1, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_red1 = (dat_new1 > test_1).float().mean()
            rej_red1 = (test_1 > cv_red1)
        if option_l in ("Gaussian", "all"):
            dat_new2 = torch.max(dat_new2, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Gaussian multiplier
            cv_norm1 = torch.kthvalue(dat_new2, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_norm1 = (dat_new2 > test_2).float().mean()
            rej_norm1 = (test_2 > cv_norm1)
        if option_l in ("Mammen", "all"):    
            dat_new3 = torch.max(dat_new3, axis = 1).values     # calculate the row-wise maximum
            # calculate the critical value with Mammen's multiplier
            cv_mm1 = torch.kthvalue(dat_new3, (N-int(N*alpha)+1)).values
            # calculate the p value 
            pval_mm1 = (dat_new3 > test_3).float().mean()
            rej_mm1 = (test_3 > cv_mm1)


    #-------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------    
    result = None
    if option_l == "Rademacher":
        result = {
            "option": "Rademacher",
            "reject": int(rej_red1.item()),
            "p_value": pval_red1.item(),
            "test_stat": test_1.item(),
            "cv": cv_red1.item(),
        }
    elif option_l == "Gaussian":
        result = {
            "option": "Gaussian",
            "reject": int(rej_norm1),
            "p_value": pval_norm1.item(),
            "test_stat": test_2.item(),
            "cv": cv_norm1.item(),
        }
    elif option_l == "Mammen":
        result = {
            "option": "Mammen",
            "reject": int(rej_mm1),
            "p_value": pval_mm1.item(),
            "test_stat": test_3.item(),
            "cv": cv_mm1.item(),
        }
    else:  # option_l == "all"
        result = {
            "Rademacher": {
                "option": "Rademacher",
                "reject": int(rej_red1.item()),
                "p_value": pval_red1.item(),
                "test_stat": test_1.item(),
                "cv": cv_red1.item(),
            },
            "Gaussian": {
                "option": "Gaussian",
                "reject": int(rej_norm1.item()),
                "p_value": pval_norm1.item(),
                "test_stat": test_2.item(),
                "cv": cv_norm1.item(),
            },
            "Mammen": {
                "option": "Mammen",
                "reject": int(rej_mm1.item()),
                "p_value": pval_mm1.item(),
                "test_stat": test_3.item(),
                "cv": cv_mm1.item(),
            },
        }
    
    return result






def gauss_fun_full(option_l,x, y, z, n, alpha,batchsize, lr,n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, drop_last1, device):    
    """
    implement the procedure without sample-splitting (see Remark 4 in Section 4.1).
    """ 
    
    # coordinatewise Gaussianization
    x = data_transform(x, x)
    y = data_transform(y, y)
    z = data_transform(z, z)

    p = x.size(1)
    q = y.size(1)
    m = z.size(1)
    
    # Concatenate response
    rep_xy = torch.cat([x, y], dim=1)
    m1=max(p,q,m) 
    thresh = (math.log(m1*n))**(1/2) * math.log(n)

    # For every column in [x, y], fit a neural network using z as predictors 
    results = [trainer("full-sample", z, rep_xy[:, i], n, batchsize,  lr, n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, device,thresh, drop_last1) for i in range(p + q)]
    
    error_all = torch.stack([res['error'] for res in results],dim=1)
    eps_all = error_all[:, 0:p]             # get residuals of x on z
    delta_all = error_all[:, p:(p+q)]       # get residuals of y on z
     
    N=5000              # number of multiplier bootstrap replications
    n_s = 0
    ####Rademacher  
    if option_l in ("Rademacher", "all"): 
        n_tilde_red = n    
        eps_t1 = eps_all[n_s:(n_s+n_tilde_red), :]
        delta_t1 = delta_all[n_s:(n_s+n_tilde_red), :] 
        #calculate test statistic 
        test_1 = torch.max(abs(eps_t1.T @ delta_t1 / np.sqrt(n_tilde_red) ))   # calculate test statistics with selected n_3^{opt} (Rademacher)
        ### generate Radamacher multipliers (N-by-n_tilde_red matrix) 
        M = generate_one_minus_one_distribution(N, n_tilde_red).float().to(device)
        
    ####Gaussian
    if option_l in ("Gaussian", "all"):
        n_tilde_norm = n     
        eps_t2 = eps_all[n_s:(n_s+n_tilde_norm), :]
        delta_t2 = delta_all[n_s:(n_s+n_tilde_norm), :]
        #calculate test statistic
        test_2 =  torch.max(abs(eps_t2.T @ delta_t2 / np.sqrt(n_tilde_norm) ))  # calculate test statistics with selected n_3^{opt} (Gaussian)
        ### generate Gaussian multipliers (N-by-n_tilde_norm matrix)
        M_norm = torch.randn(N, n_tilde_norm).to(device)
     
    ####Mammen's    
    if option_l in ("Mammen", "all"):
        n_tilde_mm = n      # choose n_3^{opt} (Mammen's) with empirical size closest to alpha  
        eps_t3 = eps_all[n_s:(n_s+n_tilde_mm), :]
        delta_t3 = delta_all[n_s:(n_s+n_tilde_mm), :]
        #calculate test statistic
        test_3 = torch.max(abs(eps_t3.T @ delta_t3 / np.sqrt(n_tilde_mm) )) # calculate test statistics with selected n_3^{opt} (Mammen's)
        ### generate Mammen's multipliers (N-by-n_tilde_mm matrix)
        M_mm = generate_mammen(N,n_tilde_mm).float().to(device)
    

    # ------------------------------
    # Calculate the critical values  
    # ------------------------------ 
    if option_l in ("Rademacher", "all"):
        n_tilde_red = torch.sqrt(torch.as_tensor(n_tilde_red)).to(device)   
    if option_l in ("Gaussian", "all"):
        n_tilde_norm = torch.sqrt(torch.as_tensor(n_tilde_norm)).to(device)  
    if option_l in ("Mammen", "all"):
        n_tilde_mm = torch.sqrt(torch.as_tensor(n_tilde_mm)).to(device)   

    dat_new1 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
    dat_new2 = torch.empty(N, p).to(device)         # generate an N-by-p matrix
    dat_new3 = torch.empty(N, p).to(device)         # generate an N-by-p matrix

    for i in range(p):
        ### Radamacher multiplier 
        if option_l in ("Rademacher", "all"): 
            dat_1 = eps_t1[:, i:i+1] * delta_t1   
            dat_1 = (dat_1 - torch.mean(dat_1, dim=0, keepdim=True)).to(device)
            dat_new1[:, i] = torch.max(abs(M @ dat_1 / n_tilde_red), axis=1)[0]  
        ### Gaussian multiplier
        if option_l in ("Gaussian", "all"):
            dat_2 = eps_t2[:, i:i+1] * delta_t2   
            dat_2 = (dat_2 - torch.mean(dat_2, dim=0, keepdim=True)).to(device)   
            dat_new2[:, i] = torch.max(abs(M_norm @ dat_2 / n_tilde_norm), axis=1)[0]   
        ### Mammen's multiplier
        if option_l in ("Mammen", "all"):
            dat_3 = eps_t3[:, i:i+1] * delta_t3   
            dat_3 = (dat_3 - torch.mean(dat_3, dim=0, keepdim=True)).to(device)   
            dat_new3[:, i] = torch.max(abs(M_mm @ dat_3 / n_tilde_mm), axis=1)[0]   
 
    if option_l in ("Rademacher", "all"):
        dat_new1 = torch.max(dat_new1, axis = 1).values     # calculate the row-wise maximum
        # calculate the critical value with Radamacher multiplier
        cv_red = torch.kthvalue(dat_new1, (N-int(N*alpha)+1)).values
        # calculate the p value 
        pval_red = (dat_new1 > test_1).float().mean()
        rej_red = (test_1 > cv_red)
    if option_l in ("Gaussian", "all"):
        dat_new2 = torch.max(dat_new2, axis = 1).values     # calculate the row-wise maximum
        # calculate the critical value with Gaussian multiplier
        cv_norm = torch.kthvalue(dat_new2, (N-int(N*alpha)+1)).values
        # calculate the p value 
        pval_norm = (dat_new2 > test_2).float().mean()
        rej_norm = (test_2 > cv_norm)
    if option_l in ("Mammen", "all"):    
        dat_new3 = torch.max(dat_new3, axis = 1).values     # calculate the row-wise maximum
        # calculate the critical value with Mammen's multiplier
        cv_mm = torch.kthvalue(dat_new3, (N-int(N*alpha)+1)).values
        # calculate the p value 
        pval_mm = (dat_new3 > test_3).float().mean()
        rej_mm = (test_3 > cv_mm)
 
  
    result = None
    if option_l == "Rademacher":
        result = {
            "option": "Rademacher",
            "reject": int(rej_red.item()),
            "p_value": pval_red.item(),
            "test_stat": test_1.item(),
            "cv": cv_red.item(),
        }
    elif option_l == "Gaussian":
        result = {
            "option": "Gaussian",
            "reject": int(rej_norm.item()),
            "p_value": pval_norm.item(),
            "test_stat": test_2.item(),
            "cv": cv_norm.item(),
        }
    elif option_l == "Mammen":
        result = {
            "option": "Mammen",
            "reject": int(rej_mm.item()),
            "p_value": pval_mm.item(),
            "test_stat": test_3.item(),
            "cv": cv_mm.item(),
        }
    else:  # option_l == "all"
        result = {
            "Rademacher": {
                "option": "Rademacher",
                "reject": int(rej_red.item()),
                "p_value": pval_red.item(),
                "test_stat": test_1.item(),
                "cv": cv_red.item(),
            },
            "Gaussian": {
                "option": "Gaussian",
                "reject": int(rej_norm.item()),
                "p_value": pval_norm.item(),
                "test_stat": test_2.item(),
                "cv": cv_norm.item(),
            },
            "Mammen": {
                "option": "Mammen",
                "reject": int(rej_mm.item()),
                "p_value": pval_mm.item(),
                "test_stat": test_3.item(),
                "cv": cv_mm.item(),
            },
        }
    
    return result



"""
Args:
    setting: (1) "S-selectn3":  splitting sample with selecting n_3^{opt}  (see Algorithm, Section 4.1).
             (2) "S-givenn3":   splitting sample but without selecting n_3^{opt} (see Remark 4, Section 4.1).
             (3) "full-sample": without sample-splitting (see Remark 4, Section 4.1).
    option_l: Rademacher / Gaussian / Mammen
    x, y, z:  data matrices with the same number of rows (sample size n) 
    alpha:    significance level in (0, 1)
Outputs:
    option:    the choose option
    reject:    1: reject the null hypothesis
    p_value:   p value
    test_stat: test statistic
    cv:        critical value 
""" 

def Cind_Gtest_py(device, setting, option_l, x, y, z, alpha, batchsize,  hidden_features1 = 128, hidden_features2 = 32, lr=0.01,n_epochs= 400, patience=30, drop_last1=True):
    n = z.size(0)
    input_features = z.size(1)
    output_features = 1 
    if setting == "S-selectn3":
        res = gauss_fun_split(setting, option_l, x, y, z, n, alpha, batchsize, lr,n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, drop_last1, device)
    if setting == "S-givenn3":
        res = gauss_fun_split(setting, option_l, x, y, z, n, alpha, batchsize, lr,n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features,drop_last1, device)
    if setting == "full-sample":
        res = gauss_fun_full(option_l, x, y, z, n, alpha, batchsize, lr,n_epochs, patience, input_features, hidden_features1, hidden_features2, output_features, drop_last1, device)
    return res  