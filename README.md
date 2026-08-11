@ EMPIRICAL REPLICATION CODE
========================================

Kyungsik Nam (ksnam@hufs.ac.kr)
Finalized: August 2026

1. Scope
--------

This folder contains the data, MATLAB code, and R code for the retained
empirical results, ordered as they appear in the papers:

  Main paper: Table 1 and Figures 1-4
  Supplement: Figures 1-5 in Section D

2. Code Folder and Software Requirements
-----------------------------------------

The default code_dir in Main_FCRGT_Code.m is:

  C:\Users\noran\Dropbox\Academic_Life(201808-Present)\2_Main_Research\0_Journal_Pub\2026_Submission\202602_Sector_CC\202509_P1_Code\Program\202608_JASA_Revision_Code

If this directory is unavailable, the script uses the location of the saved
Main_FCRGT_Code.m file. No other path needs to be edited.

The code was audited with MATLAB R2024a. It requires:

  - MATLAB
  - Statistics and Machine Learning Toolbox

Table 1 was audited with R 4.4.0. It requires the R packages fda, geigen,
and readxl.

3. Running the Code
-------------------

To reproduce Main Paper Table 1, run from this folder:

  Rscript Table1_Replication.R

It reads the LTemp sheet of Data_FCRGT.xlsx, prints the paper-ordered table,
and writes Table1_Replication_Output.csv. The expected p-values are computed from the
archived 100,000-draw Monte Carlo quantiles stored as R code in
Table1_MonteCarlo_Reference.R. The historical simulation seed was not recorded.

To reproduce all retained figures, run Main_FCRGT_Code.m. The script creates a
figures subfolder and exports all paper-linked figure components there.

To reproduce a selected figure from the MATLAB Editor:

  1. Run "Set-up and Data Loading".
  2. For Main Figures 3-4 and all Supplement Figures, also run
     "Shared Estimation for Main Figures 3-4 and Supplement Figures 1-5".
  3. Run the block bearing the desired paper figure number.

Later Supplement blocks use objects computed in Shared Estimation and, where
stated in the code, earlier paper-ordered blocks.

4. Paper-order Replication Map
------------------------------

  Main Table 1
    Script: Table1_Replication.R
    Monte Carlo reference: Table1_MonteCarlo_Reference.R
    Output: Table1_Replication_Output.csv

  Main Figure 1
    Block: Main Paper - Figure 1
    Files: GTemp_dist.png, GRP_Gr_Dist_FE2.png, GTemp_MV.png,
           GRP_Gr_FE2.png

  Main Figure 2
    Block: Main Paper - Figure 2
    Files: GW1_Shock.png, GW2_Shock.png

  Main Figure 3
    Block: Main Paper - Figure 3
    Files: SR_TR_Response_GW1.png, SR_TR_Response_GW2.png

  Main Figure 4
    Block: Main Paper - Figure 4
    Files: CC_Impact1.png, CC_Impact2.png, CC_Stat00.png

  Supplement Figure 1
    Block: Supplementary Material - Figure 1
    File:  FCRGT_Residuals.png

  Supplement Figure 2
    Block: Supplementary Material - Figure 2
    Files: FCRGT_Stationary_Scree.png,
           FCRGT_Stationary_Explained_Share.png

  Supplement Figure 3
    Block: Supplementary Material - Figure 3
    Files: FCRGT_Kappa2_GW1.png, FCRGT_Kappa2_GW2.png

  Supplement Figure 4
    Block: Supplementary Material - Figure 4
    Files: FCRGT_K_Robustness_GW1.png,
           FCRGT_K_Robustness_GW2.png

  Supplement Figure 5
    Block: Supplementary Material - Figure 5
    Files: FCRGT_K_Robustness_Climate_Impact_Densities.png,
           FCRGT_K_Robustness_Climate_Impact_Moments.png

5. Replication Data
-------------------

  Data_FCRGT.xlsx

The MATLAB and R code read this workbook directly.

6. Included R Code
------------------

  Table1_Replication.R
  Table1_MonteCarlo_Reference.R

7. Included MATLAB Functions
-----------------------------

  CC_Frac_Impacts.m
  KS_Desc_stat.m
  density_decom.m
  f_kappa_est_ci.m
  f_kappa_est_ci_GW.m
  inner_product.m
  mrsum.m

8. Regarding the tables and results in the Supplementary Appendix, see the separate file "README-SUPP" in the folder "Revision/Simulations". 
-----------------------------
