# test_cov_code

This folder contains the MATLAB code to reproduce results of **Tables S7** and **S8** in  **Section R.3** of the supplementary material.

---


## Main scripts

Running the scripts below requires sourcing the supporting functions in subfolder `mat_functions_dep/`.   


 
 
### `testcov_ind.m`
Implements the experiment for the independence test in **Section R.3** of the supplementary material.  The results are reported in **Table S7**.

### `testcov_cind.m`
Implements the experiment for the conditional independence test based on linear regressions in **Section R.3** of the supplementary material. The results are reported in **Table S8**.

 
 

## Supporting functions

### `mat_functions_dep/`
This subfolder contains MATLAB scripts that are sourced by the main scripts `testcov_ind.m` and `testcov_cind.m`.  

- `gen_gauss_XY.m`  
MATLAB implementation of the data-generation procedure for Settings S1–S4 used in the independence test experiments in **Section R.3** of the supplementary material.

 

- `generate_gauss_data.m`  
MATLAB implementation of the data-generation procedure for Setting S5 used in the independence test experiments in **Section R.3** of the supplementary material.

- `make_XY_once.m`   
MATLAB implementation of the data-generation procedure for Settings S1–S4 used in the conditional independence test experiments in **Section R.3** of the supplementary material. This script requires sourcing `gen_gauss_XY.m`.
 
- `make_XY_once_dence.m`   
MATLAB implementation of the data-generation procedure for Setting  S5 used in the conditional independence test experiments in **Section R.3** of the supplementary material.  This script requires sourcing `generate_gauss_data.m`. 

- `Ind_Gtest_mat.m`   
MATLAB implementation of the proposed independence test (Section 3 of the main paper).

- `Cind_Gtest_mat.m`   
MATLAB implementation of the proposed conditional independence test based on linear regressions (CI-Lasso, Section 4.2 of the main paper).


- `Schott_statistics.m`, `TB_statistics.m` and `Wilk_statistic.m`  
MATLAB implementations of the competing methods used in the simulation studies (Section R.3 of the supplementary material). 
 
 
### `glmnet-matlab/`
 The MATLAB code relies on the **glmnet** package, which is available at  https://hastie.su.domains/glmnet_matlab/. 
 
 
 


 