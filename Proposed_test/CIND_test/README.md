# Proposed Conditional Independence Tests 

This folder contains the code for the implementations of the proposed conditional independence tests. The details of the  proposed conditional independence tests can be found in Section 4. 



## `Cind_gaussian_fun.py`  
`Cind_Gtest_py` provides the Python implementation of the proposed conditional independence test based on nonparametric regressions (CI-FNN, Section 4.1).  

### Usage
`Cind_Gtest_py(device, setting, option_l, x, y, z, alpha, batchsize, hidden_features1 = 128, hidden_features2 = 32, lr = 0.01, n_epochs = 400, patience = 30, drop_last1 = True)`

**Inputs**  
- `device`: computational device (`"cpu"` or `"cuda"`).  
- `setting`:  the sample-splitting scheme, including
  - `"S-selectn3"`: splitting sample with selecting $n_3^{opt}$ (see Algorithm 1, Section 4.1).  
  - `"S-givenn3"`: splitting sample given $n_3$ (see Remark 4, Section 4.1).  
  - `"full-sample"`: without sample-splitting (see Remark 4, Section 4.1).  
- `option_l`: the type of multipliers (Section 5), including  
  - `"Rademacher"` uses Rademacher multipliers.  
  - `"Gaussian"` uses Gaussian multipliers.  
  - `"Mammen"` uses Mammen’s multipliers.  
  - `"all"` computes results for all three types of multipliers.  
- `x, y, z`: data matrices with the same number of rows (sample size n).  
- `alpha`: significance level in (0, 1).  
- `batchsize`: mini-batch size.  
- `hidden_features1`, `hidden_features2`: numbers of hidden units (default = 128, 32).  
- `lr`: learning rate (default = 0.01).  
- `n_epochs`: maximum number of training epochs (default = 400).  
- `patience`: early stopping patience (default = 30).  
- `drop_last1`: whether to drop the last incomplete batch (default = True).  

**Outputs**  
- `option`: the chosen `option_l`.  
- `reject`: `1` = reject the null hypothesis; `0` = do not reject.  
- `p_value`: p-value.  
- `test_stat`: test statistic.  
- `cv`: critical value.  

---

## `PCIND.R`  
`Cind_Gtest` provides the R implementation of the proposed conditional independence test based on linear regressions (CI-Lasso, Section 4.2).  

### Usage 
`Cind_Gtest(x, y, z, alpha, seed = 1, option = c("Rademacher", "Gaussian", "Mammen", "all"), N = 5000)`

**Inputs**  
- `x, y, z`: data matrices with the same number of rows (sample size n).  
- `alpha`: significance level in (0, 1).  
- `seed`: RNG seed (optional).  
- `option`: the type of multipliers (Section 5), including   
  - `"Rademacher"` uses Rademacher multipliers.  
  - `"Gaussian"` uses Gaussian multipliers.  
  - `"Mammen"` uses Mammen’s multipliers.  
  - `"all"` computes results for all three types of multipliers.  
- `N`: number of bootstrap replications (default = 5000).  

**Outputs**  
- `option`: the chosen `option`.  
- `reject`: `1` = reject the null hypothesis; `0` = do not reject.  
- `p_value`: p-value.  
- `test_sta`: test statistic.  
- `cv`: critical value.  

---

## `Cind_Gtest_mat.m`  
`Cind_Gtest_mat` provides the MATLAB implementation of the proposed conditional independence test based on linear regressions (CI-Lasso, Section 4.2).  

### Usage 
`Cind_Gtest_mat(x, y, z, alpha, option, N, seed)`

**Inputs**  
- `x, y, z`: data matrices with the same number of rows (sample size n).  
- `alpha`: significance level in (0, 1).  
- `seed`: RNG seed (optional).  
- `option`: the type of multipliers (Section 5), including    
  - `"Rademacher"` uses Rademacher multipliers.  
  - `"Gaussian"` uses Gaussian multipliers.  
  - `"Mammen"` uses Mammen’s multipliers.  
  - `"all"` computes results for all three types of multipliers.  
- `N`: number of bootstrap replications (default = 5000).  

**Outputs**  
- `option`: the chosen `option`.  
- `reject`: `1` = reject the null hypothesis; `0` = do not reject.  
- `p_value`: p-value.  
- `test_sta`: test statistic.  
- `cv`: critical value.  
