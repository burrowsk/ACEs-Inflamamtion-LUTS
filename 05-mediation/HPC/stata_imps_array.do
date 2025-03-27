/* Housekeeping */

* easier to submit the two code blocks below seperately to slurm by commenting each code block
* to submit model1 "/* Analysis: IL6 only */" comment out "/* Analysis: M=CRP; L=IL6 */" and log using "log using g-form_`exposure'_`outcome'_`imp'.log, replace"
* then submit the job to HPC using submit_imps.sh
* to submit model2 "/* Analysis: M=CRP; L=IL6 */" comment out "/* Analysis: IL6 only */" and log using "log using g-form_`exposure'_`outcome'_`imp'_model2.log, replace"
* then submit the job to HPC using submit_imps.sh

clear
set more off
set linesize 200

* Add in chkin ado path as for some reason this can't be found!
include chkin.ado, adopath

* Immputation number
args imp exposure outcome
di `imp'
di "`exposure'"
di "`outcome'"

cd /user/home/epwkb/scratch/inflam/`exposure'/`outcome'/

capture log close
log using g-form_`exposure'_`outcome'_`imp'_model2.log, replace
di "Imp number is:  " `imp'
di "Exposure is: " "`exposure'"
di "Outcome is: " "`outcome'"

use "/user/home/epwkb/scratch/inflam/stata_imps_clean.dta", clear // made using descriptives.do

keep if _mi_m ==`imp'  

/*test reg*/
*reg `outcome' `exposure'

/*
/* Analysis: IL6 only */
gformula `exposure' logIL6 `outcome' Age_blood_draw_9 BMI SDQ Constipation Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
mediation outcome(`outcome') exposure(`exposure') mediator(logIL6) ///
base_confs(Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age) ///
post_confs(Age_blood_draw_9 BMI SDQ Constipation) ///
///
commands(`outcome':logit, logIL6:regress, Age_blood_draw_9:regress, BMI:regress, SDQ:regress, Constipation:logit) ///
///
equations(`outcome': `exposure' logIL6 Age_blood_draw_9 BMI Constipation SDQ Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
///
logIL6: `exposure' Age_blood_draw_9 BMI Constipation SDQ Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
///
Age_blood_draw_9: `exposure' Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
BMI: `exposure' Age_blood_draw_9 Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
SDQ: `exposure' Age_blood_draw_9 BMI Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
Constipation: `exposure' Age_blood_draw_9 BMI SDQ Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age) ///
logOR linexp seed(69) samples(500) moreMC sim(100000)

return list
*/

/* Analysis: M=CRP; L=IL6 */
gformula `exposure' logIL6 logCRP `outcome' Age_blood_draw_9 BMI SDQ Constipation Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
mediation outcome(`outcome') exposure(`exposure') mediator(logCRP) ///
base_confs(Sex mateducation House_tenure Ethnicity Marital_status Crowding smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age) ///
post_confs(Age_blood_draw_9 BMI SDQ Constipation logIL6) ///
///
commands(`outcome':logit, logIL6:regress, logCRP:regress, Age_blood_draw_9:regress, BMI:regress, SDQ:regress, Constipation:logit) ///
///
equations(`outcome': `exposure' logCRP Age_blood_draw_9 BMI Constipation SDQ logIL6 Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
///
logCRP: `exposure' Age_blood_draw_9 BMI Constipation SDQ logIL6 Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
///
Age_blood_draw_9: `exposure' Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
BMI: `exposure' Age_blood_draw_9 Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
SDQ: `exposure' Age_blood_draw_9 BMI Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
Constipation: `exposure' Age_blood_draw_9 BMI SDQ Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age, ///
logIL6: `exposure' Age_blood_draw_9 BMI SDQ Constipation Sex i.mateducation i.House_tenure Ethnicity i.Marital_status Crowding i.smoked_pregnancy_cat matage Parity devlevel Birthweight Gestational_age) ///
logOR linexp seed(69) samples(500) moreMC sim(100000)

return list

log close