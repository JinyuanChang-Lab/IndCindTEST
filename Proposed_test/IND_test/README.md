# Proposed Independence Test

This folder contains the code for the implementations of the proposed independence test. The details of the  proposed  independence test can be found in Section 3. 

## `PIND.R`  
`Ind_Gtest` provides the R implementation of the proposed independence test.  

### Usage 
`Ind_Gtest(x, y, alpha, seed = 1, option = c("Rademacher", "Gaussian", "Mammen", "all"), N = 5000)`

**Inputs**  
- `x, y`: data matrices with the same number of rows (sample size n).  
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

## `Ind_Gtest_mat.m`  
`Ind_Gtest_mat` provides the MATLAB implementation of the proposed independence test.  

### Usage 
`Ind_Gtest_mat(x, y, alpha, option, N, seed)`

**Inputs**  
- `x, y`: data matrices with the same number of rows (sample size n).  
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
