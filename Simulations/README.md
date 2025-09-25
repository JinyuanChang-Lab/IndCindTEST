# Simulation Codes

This folder contains the code for the simulation studies reported in Section 7  of the main paper and Section R of the supplementary material.  

 
---


## IND_simulation
This subfolder contains the implementations of simulation studies for independence test in Section 7.1 (main paper) and Sections R.1-R.3 (supplement).

### IND_R

-  `mex1.R`-`mex5.R`: R scripts for the numerical studies  of the proposed independence test and the competing methods in Section 7.1 (main paper) and Section R.2 (supplement). The data generation processes follow Examples 1-5 described in Section 7.1.

- `nopower.R`: R script for the numerical studies of the proposed independence test and the competing methods  in Section R.1 (supplement).  


### IND_Matlab

- `testcov_ind.m`: MATLAB script for the numerical studies of the independence test  and the competing methods  in Section R.3 (supplement).


## CIND_simulation
This subfolder contains the implementations of simulation studies for conditional independence tests in Section 7.2 (main paper) and Sections R.3-R.4 (supplement).

### CIND_R

- `ex6.R`-`ex10.R`: R scripts for the numerical studies of the proposed conditional independence test based on linear regressions (CI-Lasso) and the competing methods in Section 7.2. The data generation processes follow Examples 6-10 described in Section 7.2.

### CIND_Matlab

- `testcov_cind.m`: MATLAB script for the numerical studies of the proposed conditional independence test based on linear regressions (CI-Lasso) and the competing methods in Section R.3 (supplement).

Note: The MATLAB code relies on the **glmnet** package, which is available at [https://hastie.su.domains/glmnet_matlab/](https://hastie.su.domains/glmnet_matlab/).

### CIND_Python

- `main_run.py`: Python script for the numerical studies of the proposed conditional independence test based on nonparametric regressions (CI-FNN) in Section 7.2. The data generation processes follow Examples 6-10 described in Section 7.2.

- `main_run_n3.py`: Python script for the numerical studies of the proposed conditional independence test based on nonparametric regressions (CI-FNN) in Section R.4 (supplement). The data generation processes follow Examples 6-10 described in Section 7.2.





