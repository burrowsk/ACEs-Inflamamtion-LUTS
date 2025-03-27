clear
set more off

use "B:\Kim\02-Inflammation-mediation\02-analysis\imp_for_stata.dta", clear // output from convert_mice_ice.R

keep _mj _mi Physical_abuse Sexual_abuse Emotional_abuse Emotional_neglect Bullying Violence_parents Substance_abuse Mental_hlth_suicide Convicted_offences Parental_separation ACE_score_classic aln qlet CRP IL6 Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume UI BMI Constipation SDQ Sex matage Parity mateducation House_tenure devlevel Ethnicity Birthweight Gestational_age Marital_status Crowding smoked_pregnancy_cat Age_blood_draw_9 daywet_7_9 bedwet_7_9 freq_urine_7_9 nocturia_7_9 urgency_7_9 void_postpone_7_9

* Recode vars ------------------------------------------------------------------

* re-value and re-label values 1/2 to 0/1 no_yes
label define no_yes 0 "No" 1 "Yes"
foreach var of varlist UI Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume Constipation {
	tab `var', nol
    recode `var' (1=0) (2=1)
    label values `var' no_yes
	tab `var', nol
}

foreach var of varlist Physical_abuse Sexual_abuse Emotional_abuse Emotional_neglect Bullying Violence_parents Substance_abuse Mental_hlth_suicide Convicted_offences Parental_separation {
	tab `var', nol
    recode `var' (1=0) (2=1)
    label values `var' no_yes
	tab `var', nol
}

label var Daywetting "Has daytime wetting"
label var Bedwetting "Has bedtime wetting"
label var UI "Has any UI"
label var Frequency "Has high frequent urination"
label var Nocturia "Has nocturia"
label var Urgency "Has urgency"
label var Voiding_postponement "Postpones voiding"
label var Voiding_volume "Low voided volume"

tab ACE_score_classic 
label var ACE_score_classic "ACEs score: Classic"
gen ACE_score = ACE_score_classic
replace ACE_score =7 if ACE_score_classic >=7 & ACE_score_classic !=.
tab ACE_score

sum CRP IL6
label var CRP "CRP @ 9 (mg/L)"
label var IL6 "IL6 @ 9 (pg/ml)"

gen logCRP = log(CRP)
gen logIL6 = log(IL6)
label var logCRP "Log CRP"
label var logIL6 "Log IL6"

tab Sex
tab Sex, nol
recode Sex (1=0) (2=1) (else=.)
label define Sex 0 "Male" 1 "Female", replace
label values Sex Sex
label var Sex "Sex (Female)"
tab Sex

sum matage
label var matage "Mothers age at delivery (yrs)"

sum Parity
label var Parity "Parity"

tab mateducation
tab mateducation, nol
recode mateducation (1=0) (2=1) (3=2)
label define mateducation_lbl 0 "A-Level or higher" 1 "O Level" 2 "Vocational or less", replace
label values mateducation mateducation_lbl
tab mateducation
label var mateducation "Maternal education: 3=vocational or less"

tab House_tenure
tab House_tenure, nol
recode House_tenure (1=0) (2=1) (else=.)
label define House_tenure_lbl 0 "Mortgaged, Owned, or Privately Rented" 1 "Rented or Other", replace
label values House_tenure House_tenure_lbl
tab House_tenure
label var House_tenure "Owned vs. renting/other housing"

sum devlevel
label var devlevel "Child developmental level 18mths (cont)"

tab Ethnicity 
tab Ethnicity, nol
recode Ethnicity (1=0) (2=1) (else=.)
label define Ethnicity_lbl 0 "White" 1 "Non-white", replace
label values Ethnicity Ethnicity_lbl
tab Ethnicity
label var Ethnicity "Ethnicity (non-white)"

sum Birthweight
label var Birthweight "Birthweight (g)"

sum Gestational_age
label var Gestational_age "Gestational age (weeks)"

tab Marital_status
tab Marital_status, nol
recode Marital_status (1=0) (2=1) (3=2)
label define marital_status_lbl 0 "Married" 1 "Single" 2 "Divorced/separated/widowed"
label values Marital_status marital_status_lbl
tab Marital_status
label var Marital_status "marital_status status"

tab Crowding
tab Crowding, nol
recode Crowding (1=0) (2=1) (else=.)
label define Crowding_lbl 0 "No Crowding" 1 "Crowding:>1"
label values Crowding Crowding_lbl
tab Crowding
label var Crowding "Crowding: >1 person per room"

tab smoked_pregnancy_cat
tab smoked_pregnancy_cat, nol
recode smoked_pregnancy_cat (1=0) (2=1) (3=2)
label define smoke 0 "None" 1 "Yes, quit early" 2 "Yes, throughout"
label values smoked_pregnancy_cat smoke
tab smoked_pregnancy_cat
label var smoked_pregnancy_cat "Smoked during pregnancy (no/quit/yes)"

tab SDQ 
sum SDQ
label var SDQ "SDQ age 8 pro-rated 0-40 score"

label var Age_blood_draw_9 "Age at F9 blood draw"
label var BMI "BMI at F@8 (kg/m2)"

label var Constipation "Has constipation"

* Order vars--------------------------------------------------------------------
order _mj _mi aln qlet ACE_score_classic ACE_score  ///
UI Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume ///
CRP IL6 logCRP logIL6 ///
Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age Age_blood_draw_9 ///
Constipation BMI SDQ


* Sorting the imputation datasets-----------------------------------------------
mi import ice, automatic clear

mi describe
mi varying

* SAVE FOR MEDIATION
drop daywet_7_9 bedwet_7_9 freq_urine_7_9 nocturia_7_9 urgency_7_9 void_postpone_7_9
save "B:\Kim\02-Inflammation-mediation\02-analysis\data\stata_imps_clean.dta", clear

* Missing data by variable type 
egen complete_cases_ACEs = rowmiss(IL6 CRP Physical_abuse Sexual_abuse Emotional_abuse Emotional_neglect Bullying Violence_parents Substance_abuse Mental_hlth_suicide Convicted_offences Parental_separation) if _mi_m==0
recode complete_cases_ACEs (0=1) (1/10=0)
label values complete_cases_ACEs cc_lbl
tab complete_cases_ACEs, missing

egen complete_cases_ls = rowmiss(IL6 CRP Age_blood_draw_9 Constipation BMI SDQ) if _mi_m==0
recode complete_cases_ls (0=1) (1/3=0)
label values complete_cases_ls cc_lbl
tab complete_cases_ls, missing

egen complete_cases_confs = rowmiss(IL6 CRP Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age) if _mi_m==0
recode complete_cases_confs (0=1) (1/33=0)
label values complete_cases_confs cc_lbl
tab complete_cases_confs, missing

egen complete_cases_out = rowmiss(IL6 CRP Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume) if _mi_m==0
recode complete_cases_out (0=1) (1/33=0)
label values complete_cases_out cc_lbl
tab complete_cases_out, missing

* Combine the cc flags to generate ultimate cc flag
egen complete_cases = rowmiss(IL6 CRP Physical_abuse Sexual_abuse Emotional_abuse Emotional_neglect Bullying Violence_parents Substance_abuse Mental_hlth_suicide Convicted_offences Parental_separation Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age Age_blood_draw_9 Constipation BMI SDQ UI Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume) if _mi_m==0
recode complete_cases (0=1) (1/33=0)
label values complete_cases cc_lbl
tab complete_cases, missing

* Cohort characteristics because doing this in R was an absolute nightmare to get the table how I wanted it to be-------
* ACES
foreach var of varlist Physical_abuse Sexual_abuse Emotional_abuse Emotional_neglect Bullying Violence_parents Substance_abuse Mental_hlth_suicide Convicted_offences Parental_separation {
mi estimate: proportion `var'
mi estimate, eform: glm `var', link(log) family(binomial) nolog
proportion `var' if complete_cases ==1
tab `var' if complete_cases ==1
}

foreach var of varlist ACE_score ACE_score_classic {
mi estimate: regress `var'
mi estimate: mean `var'
mean `var' if complete_cases ==1
sum `var' if complete_cases ==1
}

foreach var of varlist IL6 CRP {
	mi estimate: mean `var'
	mean `var' if complete_cases ==1
	sum `var' if complete_cases ==1
}

* recale vars
mi passive: gen age = Age_blood_draw_9 /12
mi passive: gen bw = Birthweight/1000

foreach var of varlist Age_blood_draw_9 age BMI SDQ matage Parity devlevel Birthweight bw Gestational_age {
	mi estimate: mean `var'
	mean `var' if complete_cases ==1
	sum `var' if complete_cases ==1
}

foreach var of varlist Constipation Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat {
	mi estimate: proportion `var'
	proportion `var' if complete_cases ==1
	tab `var' if complete_cases ==1
}

foreach var of varlist UI Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume {
	mi estimate: proportion `var'
	proportion `var' if complete_cases ==1
	tab `var' if complete_cases ==1
}

* Correlations
tetrachoric Sexual_abuse Physical_abuse Emotional_abuse Emotional_neglect Substance_abuse Mental_hlth_suicide Violence_parents Parental_separation Bullying Convicted_offences if _mi_m ==0, stats(rho p) star(0.05)
* broadly consistent with some fluctuations  - there are missing ACEs in ARs table

tetrachoric UI Daywetting Bedwetting Frequency Nocturia Urgency Voiding_postponement Voiding_volume if _mi_m ==0, stats(rho p) star(0.05)

/*
* Check cell counts for the low prevalence ACEs (doesn't matter anymore as I am no longer exploring individual ACEs)
mi xeq: tab Sexual_abuse UI // okay
mi xeq: tab Sexual_abuse Daywetting // cells =0
mi xeq: tab Sexual_abuse Bedwetting //  cells =0
mi xeq: tab Sexual_abuse Frequency // okay
mi xeq: tab Sexual_abuse Nocturia // okay
mi xeq: tab Sexual_abuse Urgency //  cells =0
mi xeq: tab Sexual_abuse Voiding_postponement // okay
mi xeq: tab Sexual_abuse Voiding_volume //  cells =0

mi xeq: tab Bullying UI // okay
mi xeq: tab Bullying Daywetting // cells =0
mi xeq: tab Bullying Bedwetting //  cells =0
mi xeq: tab Bullying Frequency // okay
mi xeq: tab Bullying Nocturia // okay
mi xeq: tab Bullying Urgency //  cells =0
mi xeq: tab Bullying Voiding_postponement // okay
mi xeq: tab Bullying Voiding_volume //  cells =0

* Test a regression
mi estimate, noisily or: logistic Voiding_volume Sexual_abuse Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age

/*
Error message:
mi estimate: omitted terms vary
    The set of omitted variables or categories is not consistent between m=1 and m=31; this is not allowed.  To identify varying sets, you can use mi xeq to
    run the command on individual imputations or you can reissue the command with mi estimate, noisily
*/
*/

* Testing regressions
mi estimate, or: logistic UI ACE_score
mi estimate, or: logistic UI ACE_score Sex 
mi estimate, or: logistic UI ACE_score Sex i.mateducation 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight 
mi estimate, or: logistic UI ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age

mi estimate, eform: regress logIL6 ACE_score
mi estimate, eform: regress logIL6 ACE_score Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age

mi estimate, or: logistic UI logIL6
mi estimate, or: logistic UI logIL6 ACE_score Age_blood_draw_9 BMI SDQ Constipation Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age

mi estimate, or: logistic UI logCRP
mi estimate, or: logistic UI logCRP logIL6 Physical_abuse Age_blood_draw_9 BMI SDQ Constipation Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age

* END