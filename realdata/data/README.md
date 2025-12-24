# data
This folder contains the raw data and the processed data used in the real data analysis in **Section Q** of the supplementary material.

## `raw_data/`


This subfolder contains the raw daily stock price data obtained from the WRDS website https://wrds-www.wharton.upenn.edu/ by navigating to:

Get Data $\to$ CRSP $\to$ Stock / Security Files $\to$ Daily Stock File.

On the Daily Stock File page https://wrds-www.wharton.upenn.edu/pages/get-data/center-research-security-prices-crsp/annual-update/stock-security-files/daily-stock-file/, the following steps are performed:  
***Step 1***  Set the sample period from 2016-01-01 to 2018-12-31 for the before COVID-19 period.  
***Step 2***  Under the Autocomplete option, select “ticker” to specify the stocks. Then click “Browse…” to upload the plain-text ticker list  `<sector>.txt` (stored in the `raw_data/classification` folder) for each industrial sector.  
***Step 3*** Select PERMNO, TICKER, and PRC as the output variables, defined as follows:  
    - PERMNO: the permanent security identifier assigned by CRSP, used as the unique identifier in our empirical analysis;  
    - TICKER: the stock ticker name (as specified in the ticker list `<sector>.txt`);   
    - PRC: the stock price (either the transaction price or the bid/ask average).  
***Step 4*** Select “.csv” as the output format, and submit the query to retrieve the data. The output file is named `68<sector>.csv` (e.g., `68CD.csv` ) for each sector. 

To retrieve data for the during/after COVID-19 period, set the sample period in ***Step 1*** from 2020-01-01 to 2022-12-31 and follow the same steps. These data are saved in files named `02<sector>.csv` (e.g., `02CD.csv`) for each sector.


**Table 1. Raw data files**

| Sector name | 2016/01/01–2018/12/31 (before COVID-19) | 2020/01/01–2022/12/31 (during/after COVID-19) |
|:---|:---:|:---:|
| Consumer Discretionary (CD) | `68CD.csv` | `02CD.csv` |
| Communication Services (CS) | `68CS.csv` | `02CS.csv` |
| Consumer Staples (CSt) | `68CSt.csv` | `02CSt.csv` |
| Energy (Eng) | `68Eng.csv` | `02Eng.csv` |
| Financials (Fin) | `68Fin.csv` | `02Fin.csv` |
| Health Care (HC) | `68HC.csv` | `02HC.csv` |
| Industrials (Ind) | `68Ind.csv` | `02Ind.csv` |
| Information Technology (IT) | `68IT.csv` | `02IT.csv` |
| Materials (Mat) | `68Mat.csv` | `02Mat.csv` |
| Real Estate (RE) | `68RE.csv` | `02RE.csv` |
| Utilities (Uti) | `68Uti.csv` | `02Uti.csv` |

 
 

##  `processed_data/`

 
This subfolder contains the daily stock returns after data preprocessing by using `data_processing.R`. 

The data preprocessing procedures include the following steps:  
(i) aligning the time series of different stocks by trading dates;   
(ii) removing stocks with missing observations;  
(iii) computing log returns;  
(iv) testing for stationarity with the ADF unit-root test. 


The outputs are saved as:

**Table 2. Processed data files**

| Sector name | 2016/01/01–2018/12/31 (before COVID-19) | 2020/01/01–2022/12/31 (during/after COVID-19) |
|:---|:---:|:---:|
| Consumer Discretionary (CD) | `68return_CD.csv` | `02return_CD.csv` |
| Communication Services (CS) | `68return_CS.csv` | `02return_CS.csv` |
| Consumer Staples (CSt) | `68return_CSt.csv` | `02return_CSt.csv` |
| Energy (Eng) | `68return_Eng.csv` | `02return_Eng.csv` |
| Financials (Fin) | `68return_Fin.csv` | `02return_Fin.csv` |
| Health Care (HC) | `68return_HC.csv` | `02return_HC.csv` |
| Industrials (Ind) | `68return_Ind.csv` | `02return_Ind.csv` |
| Information Technology (IT) | `68return_IT.csv` | `02return_IT.csv` |
| Materials (Mat) | `68return_Mat.csv` | `02return_Mat.csv` |
| Real Estate (RE) | `68return_RE.csv` | `02return_RE.csv` |
| Utilities (Uti) | `68return_Uti.csv` | `02return_Uti.csv` |
