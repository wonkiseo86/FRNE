@ Empirical Program Code Descriptions for “Functional Regression with Persistent Nonstationarity and Error Contamination”
                                                                                                                                                   September 2025
                                                                                                                                                      Kyungsik Nam

1. Main Analysis
The primary empirical analysis is implemented in the file Main_FCRGT_Code.m, which reproduces the results presented in Section 5 (Figures 1–4) of the paper.

2. Data Files
TPDN_g.mat: Contains the density-valued functional time series of land temperature anomalies from January 1951 to December 2019.
GRP_den2_full.mat: Contains the density-valued functional time series of temperature-related GRP growth rates over the period 1951–2019.

3. Variance Ratio (VR) Test (Table 1)
To replicate the VR test reported in Table 1, use the R script FCRGT_Nonstat_Test.r together with the Excel data file Data_FCRGT.xlsx.

4. Comp_pvals: Code for computing p-values for the VR test statistics.

5. Comp_sim: Code used for simulation experiments. 

6. Contact Information
For questions regarding the program code or data, please contact: Kyungsik Nam (ksnam@hufs.ac.kr) or Won-Ki Seo (won-ki.seo@sydney.edu.au).
