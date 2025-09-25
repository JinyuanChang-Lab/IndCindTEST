# Real Data Analysis 
We investigate the dependence and conditional dependence structures in the S&P 500 stocks in the real data analysis (Section Q in the supplementary material). 


## Data
 We collect the daily stock prices (485 stocks) for 11 industrial sectors over two non-overlapping windows: 
  
  (a)	from 1 January 2016 to 31 December 2018 (before COVID-19 period);
  
  (b)	from 1 January 2020 to 31 December 2022 (during/after COVID-19 period). 

The sector classification follows the Global Industry Classification Standard, which includes Consumer Discretionary (CD), Communication Services (CS), Consumer Staples (CSt), Energy (Eng), Financials (Fin), Health Care (HC), Industrials (Ind), Information Technology (IT), Materials (Mat), Real Estate (RE), and Utilities (Uti). The datasets are stored in two subfolders: "2016-2018data" and "2020-2022data."


- **2016-2018data**: daily stock price data from 1 January 2016 to 31 December 2018 (before COVID-19 period), stored in files named "68\<sector\>.csv" (e.g., 68CD.csv). 

- **2020-2022data**: daily stock price data from 1 January 2020 to 31 December 2022 (during/after COVID-19 period), stored in files named "02\<sector\>.csv" (e.g., 02CD.csv).


Data preprocessing for each sector in both periods was conducted using the R script "data_processing.R".  The daily stock returns can be obtained after preprocessing, which are saved as "68\<sector\>_return.csv" (e.g., 68CD_return.csv) for the before COVID-19 period and "02\<sector\>_return.csv" (e.g., 02CD_return.csv) for the during/after COVID-19 period. These files are available in the subfolder **realdata_code**.  


## Code
The subfolder **realdata_code** contains the R and Python code for the real data analysis.  To run the code, the files 68\<sector\>_return.csv and 02\<sector\>_return.csv must be placed in the same directory.  

- `ind-realdat_pval.R`: R script for computing p-values of the  proposed independence test on the real data.

- `cind-realdat_cpval.R`: R script for computing p-values of the  proposed conditional independence test based on linear regressions (CI-Lasso) and the competing methods on the real data.

- `realdata_nn.ipynb`: Python script for computing p-values of the  proposed conditional independence test based on nonparametric regressions (CI-FNN) on the real data.
