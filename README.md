# Testing Independence and Conditional Independence in High Dimensions via  Coordinatewise Gaussianization

## Introduction

We propose new statistical tests, in high-dimensional settings, for testing the independence of two random vectors and their conditional independence given a third
random vector. The key idea is simple, i.e., we first transform each component variable to the standard normal via its marginal empirical distribution, and we then test for independence and conditional independence of the transformed random vectors using appropriate $L_{\infty}$-type test statistics. While we are testing some necessary conditions of the independence or the conditional independence, the new tests outperform the 13 frequently used testing methods in a large scale simulation comparison. The advantage of the new tests can be summarized as follows: (i) they do not require any moment
conditions, (ii) they allow arbitrary dependence structures of the components among the random vectors, and (iii) they allow the dimensions of random vectors to diverge at the exponential rates of the sample size. The critical values of the proposed tests are determined by a computationally efficient multiplier bootstrap procedure. Theoretical
analysis shows that the sizes of the proposed tests can be well controlled by the nominal significance level, and the proposed tests are also consistent under certain local alternatives. The finite sample performance of the new tests is illustrated via extensive simulation studies and a real data application.


---

## Repository Structure

### **Proposed_test**  
This folder contains the code for the implementations of the proposed independence and conditional independence tests.  
- **IND_test**: implementation of the proposed independence test (Section 3).  
- **CIND_test**: implementation of the proposed conditional independence tests (Section 4).  

### **Simulations**  
This folder contains the code for simulation studies in Section 7 of the main paper and Section R of the supplementary material.  
- **IND_simulation**: implementation of simulation studies for independence tests.  
  - **IND_R**: R code for  the numerical studies in Section 7.1 (main paper) and Sections R.1-R.2 (supplement).  
  - **IND_Matlab**: MATLAB code for  the numerical studies in Section R.3 (supplement).  

- **CIND_simulation**: implementation of simulation studies for conditional independence tests.  
  - **CIND_R**: R code for the numerical studies (CI-Lasso and the competing methods) in Section 7.2.
  - **CIND_Python**: Python code for the numerical studies (CI-FNN) in Section 7.2 and Section R.4 (supplement).  
  - **CIND_Matlab**: MATLAB code for the numerical studies (CI-Lasso and the competing methods) in Section R.3 (supplement).  

### **Real_data**  
 This folder contains the data and code for the real data analysis in Section Q of the supplementary material.   The data are publicly available in the Wharton Research Data Services (WRDS) database at https://wrds-www.wharton.upenn.edu/.


- **2016-2018data**: daily stock price data from 1 January 2016 to 31 December 2018 (before COVID-19 period). 
- **2020-2022data**: daily stock price data from 1 January 2020 to 31 December 2022 (during/after COVID-19 period).

- **realdata_code**: R and Python code used to conduct the real data analysis. 
