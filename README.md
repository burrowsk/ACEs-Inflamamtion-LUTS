# ACEs-Inflammation-LUTS
Code to accompany the paper: Adverse childhood experiences and lower urinary tract symptoms in adolescence: the mediating role of inflammation

## Order of analysis

### 1. Make Phenotype File
contact kimberley.burrows[@bristol.ac.uk] for advice on creating the initial phenotype file from ALSPAC syntax scripts

### 02-make-imputations
Run on local machine (RStudio):  
- 01-Script_1.R
- 02_Script_2.R

[These two scripts format, clean, and makes the ACE constructs]

Run using HPC:  
- submit_imputations.sh [submission file to submit 03_Script_3_HPC.R to Slurm]
- 03_Script_3_HPC.R [Formats, imputes data, checks imputations]
This script sources the following scripts:
- 03a_my_vars.R
- 03b_vars_to_drop.R
- 03c_names.R
- 03d_logreg_2.R
- 03e_prop_plot_bin.R
- 03f_prop_plot_cat.R

[These are scripts that contain various vectors or functions used to complete *03_Script_3_HPC.R*]

### 03-descriptives
Run using R:
- convert_mice_ice.R [converts the mice format imputed data to ICE format (easier to read into Stata)

Run using Stata:
- descriptives_stata.do [descriptive statistics for Table S5]

### 04-regressions
Run using local machine (RStudio):
- regressions.R [regressions for Tables 1, 2, and S8 using imputed data]
- regressions_cc.R [regressions for Tables S9 and S10 using complete case data]

### 05-mediation
Run using HPC:
- submit_submit.sh
- submit_imps.sh
[these two scripts are used to submit the multiple mediations models (1 exposure * 8 outcomes * 80 imputed datasets * 2 models)
- stata_imps_array.do [this is the syntax to perform the mediation analysis for model 1 and 2]
- search_mem_issues.sh [search the Stata log files for error messages due to memory issues]
- app_logs.sh [append the 80 log files into 1 file for each of the exposure-outcome models]

### 05-mediation/complete_case
Run using HPC:
- submit_submit.sh
- submit_imps.sh
[these two scripts are used to submit the multiple mediations models (1 exposure * 8 outcomes * 2 models)
- stata_imps_array.do [this is the syntax to perform the mediation analysis for model 1 and 2]
- search_mem_issues.sh [search the Stata log files for error messages due to memory issues]
- app_logs.sh [append the 80 log files into 1 file for each of the exposure-outcome models]
- change_filenames.sh [can't quite remember why this was needed]

### 05-mediation
- rrules_ace_score_luts.R
- rrules_ace_score_luts_mod2.R

[scripts to perform Rubin's Rules for model 1 and 2 of imputed data and create Table 3 and S11]
