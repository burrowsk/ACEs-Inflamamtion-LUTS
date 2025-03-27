# Mediation of inflammation in the associations between ACEs and LUTS
#===============================================================================
# This is an adaptation of Script_3.R created by Lotte Houtepen
# Kimberley Burrows
#===============================================================================

#-------------------------------------------------------------------------------
# Admin ------------------------------------------------------------------------

rm(list = ls())
#renv::status()
#renv::snapshot()

# Library
library(foreign)
library(haven)
library(readxl)
library(matrixStats)
library(tableone)
library(mice)
library(dplyr)
library(VIM)
library(purrr)

# Notes: -----------------------------------------------------------------------
# Uses data output from Script1.R then Script2.R
# Requires: SI_data1_ACE_definitions_overtime.xlsx alspacKids_ACE_0_8_imputation.RData
# Subset to sample with data on IL6 and CRP (N=4745)
#-------------------------------------------------------------------------------


#===============================================================================
# Options
#===============================================================================

### 0. Location files
loc_inp = '~/scratch/aces/data/'
loc_out = '~/scratch/aces/data/' #  output data
setwd(loc_out)

### 1. Time-period for ACEs
timeperiod <- c(0, 8)
fileid_out <-
  paste0('alspacKids_ACE_', timeperiod[1], '_', timeperiod[2])

### 2. load data
# a. excel file with description variables
adv_description <-
  data.frame(readxl::read_excel(
    paste0(loc_inp, 'SI_data1_ACE_definitions_overtime.xlsx'),
    sheet = 'ACE variables',
    col_names = TRUE
  ))
## b. output from Script_2.R
load(paste0(loc_inp, fileid_out, "_imputation.RData"))

## Log file initiation for output
catf <-
  function(...,
           file = paste0(loc_out, fileid_out, '_imputation_log.txt'),
           append = TRUE) {
    cat(..., file = file, append = append)
  }


#===============================================================================
# Code block 1: Add in my variables to the ACEs
# I have confounders that are ACE aux variables; I want to remove these from the ACEs
#===============================================================================

source("../scripts/02-ACEs/03a_my_vars.R") # vector of column names for non-ACE specific data

## Read in non-ACE variables (from Stata)
other_vars = read_dta(paste0(loc_inp, 'pheno_other_vars.dta'))
other_vars <- zap_missing(other_vars)
other_vars <- zap_formats(other_vars)
other_vars <- zap_label(other_vars)

## Convert bin/cat to factors
other_vars = other_vars %>% mutate(across(
  c(
    "sex",
    "mateducation",
    "marital",
    "smoked_pregnancy_cat",
    "daywet_14",
    "bedwet_14",
    "freq_urine_14",
    "nocturia_14",
    "urgency_14",
    "void_postpone_14",
    "void_volume_14",
    "ui_14",
    "constipation_8",
    "housetenure",
    "ethnicity",
    "crowding",
    "daywet_7_9",
    "bedwet_7_9",
    "freq_urine_7_9",
    "nocturia_7_9",
    "urgency_7_9",
    "void_postpone_7_9",
    "signs_urine_7_9",
    "day_soiling_7_9",
    "constipation_7"
  ),
  as_factor
))

other_vars <-
  zap_labels(other_vars) # zap labels for continuous vars

## Merge ACEs with non-ACE variables
names.use <- names(other_vars)[(names(other_vars) %in% my_vars)]
other_vars <- other_vars[, names.use]
str(other_vars, list.len = ncol(other_vars)) # Check class of variables
other_vars <- droplevels(other_vars)
alspacKids_ACE_data_2018 <-
  merge(
    alspacKids_ACE_data_2018,
    other_vars,
    by = c("aln", "qlet"),
    all.y = T
  )

### Drop variables in ACEs data that are not required
## [variables that are confounders]
## [ACE categories - no longer required]
## [glyca]

source("../scripts/02-ACEs/03b_vars_to_drop.R") # vector of column names to drop from ACE dataset
alspacKids_ACE_data_2018 = alspacKids_ACE_data_2018[, !(names(alspacKids_ACE_data_2018) %in% drop_vars)]


#===============================================================================
# Code block 2: Subset to those with inflammation data and describe
# N should = 4,745
#===============================================================================

alspacKids_ACE_data_2018 <-
  alspacKids_ACE_data_2018 %>% filter_at(vars(il6, crp), all_vars(!is.na(.)))
alspacKids_ACE_data_2018 <-
  droplevels(alspacKids_ACE_data_2018) # 4,745 individuals

str(alspacKids_ACE_data_2018,
    list.len = ncol(alspacKids_ACE_data_2018))
#sapply(alspacKids_ACE_data_2018, class)

## Describe
other_vars %>% map(summary)

## make log il6 and crp
summary(alspacKids_ACE_data_2018$il6)
summary(alspacKids_ACE_data_2018$crp)
alspacKids_ACE_data_2018$il6 <- log(alspacKids_ACE_data_2018$il6)
alspacKids_ACE_data_2018$crp <- log(alspacKids_ACE_data_2018$crp)
summary(alspacKids_ACE_data_2018$il6)
summary(alspacKids_ACE_data_2018$crp)

#===============================================================================
# Code block 3: Split gender
#===============================================================================

catf('\n\n', paste0(Sys.time()),
     'Code block 3: Split gender\n')

alspacBoys_ACE_data_2018 <-
  alspacKids_ACE_data_2018[grepl('Male', alspacKids_ACE_data_2018[['sex']]), ]

alspacGirls_ACE_data_2018 <-
  alspacKids_ACE_data_2018[grepl('Female', alspacKids_ACE_data_2018$sex), ]


#===============================================================================
# Code block 4: Determine variables for imputation
# For imputation to work, auxiliary ACE variables need to have at least 50 people in all levels.
# This needs to be true for boys and girls as these are imputed separately
#===============================================================================

catf('\n\n',
     paste0(Sys.time()),
     'Code block 2: Determine variables for imputation\n')

imp_standard <-
  grep('org', names(alspacKids_ACE_data_2018), value = T)

## check for factor variables:
table(!sapply(alspacKids_ACE_data_2018[, imp_standard], function(x)
  is.numeric(x) | is.integer(x)))

imp_factor <-
  imp_standard[!sapply(alspacKids_ACE_data_2018[, imp_standard], function(x)
    is.numeric(x) | is.integer(x))]

catf(
  '\nThere are',
  length(imp_standard),
  'imputation variables,',
  'of which',
  length(imp_factor),
  'are categorical.\n'
)

### i. identify ACE imputation variables with less than 50 people in any factor level:
imp_ACE <-
  imp_factor[gsub('_org', '', imp_factor) %in% adv_description$variable_unique]

imp_ACE_less_n50 <-
  imp_ACE[sapply(alspacBoys_ACE_data_2018[, imp_ACE], function(x)
    any(table(factor(x)) < 50)) |
      sapply(alspacGirls_ACE_data_2018[, imp_ACE], function(x)
        any(table(factor(x)) < 50))]
catf(
  '\n',
  length(imp_ACE_less_n50),
  'of the',
  length(imp_ACE),
  'categorical ACE variables,',
  'have <50 people in all factor levels in either boys or girls.\n'
)

### ii. Identify which ACE imputation variables, do fullfill criteria for binary version
imp_ACE_use_binary <-
  imp_ACE_less_n50[colSums(alspacBoys_ACE_data_2018[, gsub('_org', '', imp_ACE_less_n50)], na.rm =
                             T) >= 50 &
                     colSums(alspacGirls_ACE_data_2018[, gsub('_org', '', imp_ACE_less_n50)], na.rm =
                               T) >= 50]

### iii. If applicable remove imp_ACE_less_n50 in both original and binary version
rm_imp <-
  imp_ACE_less_n50[!imp_ACE_less_n50 %in% imp_ACE_use_binary]
length(imp_standard)

catf(length(rm_imp),
     'of these',
     length(imp_ACE_less_n50),
     'variables need to be removed')

if (length(rm_imp) > 0) {
  imp_standard <- imp_standard[!imp_standard %in% rm_imp]
}
length(imp_standard)

### iv. If applicable, make sure binary version is used imp_ACE_use_binary
sum(grepl('_org', imp_standard))

catf(
  ', but',
  length(imp_ACE_use_binary),
  'variables can be included in their dichotomised version:\n',
  imp_ACE_use_binary,
  '\n'
)

if (length(imp_ACE_use_binary) > 0) {
  imp_standard[imp_standard %in% imp_ACE_use_binary] <-
    gsub('_org', '', imp_standard[imp_standard %in% imp_ACE_use_binary])
}
sum(grepl('_org', imp_standard))

catf(
  '\nSo the total number of imputation variables is',
  length(imp_standard),
  'of which',
  sum(grepl('_org', imp_standard)),
  'are included in their original non-dichotomised format.\n'
)

## clean up
rm(
  imp_ACE,
  imp_ACE_less_n50,
  imp_ACE_use_binary,
  imp_factor,
  rm_imp,
  names.use,
  other_vars,
  adv_description
)


#===============================================================================
# Code block 5: Check variables and explore missingness
#===============================================================================

catf('\n\n',
     paste0(Sys.time()),
     'Variable check and explore missing data\n')

## Define all substantive model variables:
model_vars <- c(
  "ACEscore_classic_0_8yrs",
  "physical_abuse_0_8yrs",
  "sexual_abuse_0_8yrs",
  "emotional_abuse_0_8yrs",
  "emotional_neglect_0_8yrs",
  "bullying_0_8yrs",
  "violence_between_parents_0_8yrs",
  "substance_household_0_8yrs",
  "mental_health_problems_or_suicide_0_8yrs",
  "parent_convicted_offence_0_8yrs",
  "parental_separation_0_8yrs",
  "crp",
  "il6",
  "ui_14",
  "daywet_14",
  "bedwet_14",
  "freq_urine_14",
  "nocturia_14",
  "urgency_14",
  "void_postpone_14",
  "void_volume_14",
  "bmi_8",
  "constipation_8",
  "sdq_8",
  "sex",
  "matage",
  "parity",
  "mateducation",
  "housetenure",
  "devlevel",
  "ethnicity",
  "birthweight",
  "gest_age",
  "marital",
  "crowding",
  "smoked_pregnancy_cat",
  "age_bloods_9"
)

## Define example model variables:
model_vars_light <- c("ACEscore_classic_0_8yrs",
                      "bmi_8",
                      "constipation_8",
                      "sdq_8",
                      "il6",
                      "ui_14")

## Patterns of missing data for model_vars light: too many otherwise to plot
pdf(paste0(loc_out, '01-alspacKids-missing.pdf'))
md.pattern(alspacKids_ACE_data_2018[, model_vars_light], rotate.names = TRUE)
dev.off()

## Patterns for all model variables for output to excel
pattern <-
  md.pattern(alspacKids_ACE_data_2018[, model_vars], plot = FALSE)
pattern <- as.data.frame(pattern)
write.csv(pattern,
          file = paste0(loc_out, "01-patterns.csv"),
          row.names = T)

## Another way to see missing data counts
misscount <- numeric(nrow(alspacKids_ACE_data_2018))
for (i in 1:nrow(alspacKids_ACE_data_2018)) {
  misscount[i] <-
    VIM::countNA(alspacKids_ACE_data_2018[i, model_vars])
}

## Display the counts in a table
table(misscount)
round(table(misscount) / sum(table(misscount)) * 100, 2)

catf(
  '\n\n',
  paste0(Sys.time()),
  'Explore missing data: counts\n',
  table(misscount),
  "\n\n",
  round(table(misscount) / sum(table(misscount)) * 100, 2),
  "/n/n"
)

rm(model_vars,
   model_vars_light,
   misscount,
   i,
   pattern)


#===============================================================================
# Code block 6: Subset data for imputation
# 1. Determine all ACE variables and set in correct order
# 2. Check required variables are in dataset
# 3. Change variable names
#===============================================================================

### 1.
## Vector of all ACE measures needed/available
ACEmeasures <-
  c(
    "ACEscore_classic",
    #"ACEcat_classic",
    "physical_abuse",
    "sexual_abuse",
    "emotional_abuse",
    "emotional_neglect",
    "bullying",
    "violence_between_parents",
    "substance_household",
    "mental_health_problems_or_suicide",
    "parent_convicted_offence",
    "parental_separation",
    "ACEscore_extended",
    #"ACEcat_extended",
    "social_class",
    "financial_difficulties",
    "neighbourhood",
    #"social_support_child",
    "social_support_parent",
    #"violence_between_child_and_partner",
    "physical_illness_child",
    "physical_illness_parent",
    "parent_child_bond"
  ) # ACEs in correct order / ACEs commented out as not required or not enough data

ACEs <- ACEmeasures[c(2:11, 13:19)] # 21 for all aces
classicACEs <- ACEs[1:10]

ACEmeasures <-
  grep(paste0(ACEmeasures, collapse = '|'),
       names(alspacKids_ACE_data_2018),
       value = T)
ACEs <-
  grep(paste0(ACEs, collapse = '|'),
       names(alspacKids_ACE_data_2018),
       value = T)
classicACEs <-
  grep(paste0(classicACEs, collapse = '|'),
       names(alspacKids_ACE_data_2018),
       value = T)

catf('\n\n',
     paste0(Sys.time()),
     'Code block 6: Create predictor matrix\n')

### 2. Check variables are in file
all(imp_standard %in% names(alspacKids_ACE_data_2018))
all(ACEmeasures %in% names(alspacKids_ACE_data_2018))
all(my_vars %in% names(alspacKids_ACE_data_2018))
my_vars [!my_vars %in% names(alspacKids_ACE_data_2018)]

table(unique(c(imp_standard, ACEmeasures, my_vars)) %in% names(alspacKids_ACE_data_2018)) # all 50 ACE vars present; all 80 vars present

### 3. Subset data and variable name change
dat <-
  alspacKids_ACE_data_2018[, unique(c(ACEmeasures, imp_standard, my_vars))]
dat_boys <-
  alspacBoys_ACE_data_2018[, unique(c(ACEmeasures, imp_standard, my_vars))]
dat_girls <-
  alspacGirls_ACE_data_2018[, unique(c(ACEmeasures, imp_standard, my_vars))]
source("../scripts/02-ACEs/03c_names.R")

## Check data (names and class)
str(dat,
    list.len = ncol(dat))
#sapply(dat, class)

#-------------------------------------------------------------------------------
### There is a stray variable!
# t5412_org is a constant (no variation) and remains even after removing vars <50, probably because there is only 1 cat (>50)
# drop this!
dat <- dat[, !(names(dat) %in% c("t5412_org"))]
dat_boys <- dat_boys[, !(names(dat_boys) %in% c("t5412_org"))]
dat_girls <- dat_girls[, !(names(dat_girls) %in% c("t5412_org"))]
#-------------------------------------------------------------------------------

## Write complete case data
write.dta(dat,
          file = paste0(loc_out, "02-study-sample-data.dta"))
saveRDS(dat,
        file = paste0(loc_out, "02-study-sample-data.RData"))


#===============================================================================
# Code block 7: Create prediction matrix
# 1. Initialise MI data to create prediction matrix and methods
# 2. Adapt pred_mat
# 3. Adapt meth
#===============================================================================

## User amended logreg function (MICE) to increase n iterations n=in glm.fit
#source("../scripts/02-ACEs/03d_logreg_2.R") # n iterations = 600

### 1. Initialise the MI data -------------------------------------------------------
ini <-
  mice(
    dat,
    maxit = 0,
    pri = T,
    defaultMethod = c("pmm", "logreg", "polyreg", "polr")
  )
pred_mat <- ini$pred
meth <- ini$meth
length(meth)
dim(pred_mat)
# LoggedEvent = qlet constant  - this is not a problem it is not imputed

### 2. Adapt the prediction matrix and methods
# i. ACEs will impute ACEs + ACE auxs but not my_vars
# ii. ACE score will impute my_vars
# iii. ACE aux will impute ACE auxs + ACEs + my_vars
# iv. my_vars will impute ACEs + ACE auxs + my_vars
# v. Mat weight and BMI will not mutually impute each other
# vi. UI is passively imputed
# vii. UI will not impute other stuff
# viii. ACE scores are passively imputed
# ix. Sex is not imputed nor used to impute
# x. aln + qlet needed in datset but not imputation

## Get column numbers [https://stackoverflow.com/questions/20369145/refer-to-range-of-columns-by-name-in-r]
coln <- function(X){
  y <- rbind(seq(1,ncol(X)))
  colnames(y) <- colnames(X)
  rownames(y) <- "col.number"
  return(y)}
select
coln(pred_mat)

## Set to 0 variables that do not or are not imputed e.g passive imputation
pred_mat[, "Sex"] <- 0
pred_mat["Sex",] <- 0
pred_mat[, grep('^ACE', colnames(pred_mat))] <- 0
pred_mat[grep('^ACE', rownames(pred_mat)),] <- 0
pred_mat[, "aln"] <- 0
pred_mat["aln",] <- 0
pred_mat[, "qlet"] <- 0
pred_mat["qlet",] <- 0
pred_mat["Maternal_weight","Maternal_BMI"] <-0
pred_mat["Maternal_BMI","Maternal_weight"] <-0
pred_mat[68:75,68:75] <-0 # LUTS @ 7_9 will not impute LUTS @ 7_9
pred_mat[44:50,44:50] <-0 # LUTS @ 14 will not impute LUTS @ 14
pred_mat[31:37,31:37] <-0 # Partner of child domestic violence won't impute each other
pred_mat[, "UI"] <- 0
pred_mat["UI",] <- 0
pred_mat["CRP",] <- 0
pred_mat["IL6",] <- 0

### Deprecated: 2. Adapt the prediction matrix in excel
#-------------------------------------------------------------------------------
## Manually adapt the pred_mat in excel (as many vars)
write.csv(pred_mat , file = "03-pred-matrix.csv", row.names = T)
## Load pred_mat once adapted
#### pred_mat_org <- pred_mat
#### pred_mat <- read.csv("pred_matrix/pred_matrix_luts_b.csv")
#### rownames(pred_mat) <- pred_mat$X
#### pred_mat$X <- NULL
#### pred_mat = apply(pred_mat, 2, function(x)
####   setNames(as.numeric(paste(x)), names(x))) # all numeric (not integers)
## for troubleshooting
#### pred_mat_org <- as.matrix(pred_mat)
#### pred_mat <- read.csv("pred_matrix.csv")
#### identical(pred_mat, pred_mat_org)
#### all.equal(pred_mat, pred_mat_org)
#### dimnames(pred_mat)
#### dimnames(pred_mat_org)

## No longer needed as coded above
#-------------------------------------------------------------------------------

### 2. Adapt the methods

## rename for methods
ACEmeasures<-#ACEs in correct order
  c("ACE_score_classic", #"ACEcat_classic",
    "Physical_abuse", "Sexual_abuse", "Emotional_abuse", "Emotional_neglect",
    "Bullying","Violence_parents", "Substance_abuse",
    "Mental_hlth_suicide","Convicted_offences",
    "Parental_separation",
    "ACE_score_extended", #"ACEcat_extended",
    "Social_class", "Financial_difficulties", "Neighbourhood",
    #"social_support_child",
    "Social_support_parents",
    #"violence_between_child_and_partner",
    "Physical_illness_child",
    "Physical_illness_parents","Parent_child_bond")
ACEs<-ACEmeasures[c(2:11,13:19)] # 19 for all aces
classicACEs<-ACEs[1:10] # ACEs in correct order

ini <-
  mice(
    dat,
    maxit = 0,
    pri = T,
    predictorMatrix = pred_mat,
    defaultMethod = c("pmm", "logreg", "polyreg", "polr")
  )
meth <- ini$meth # formulas to be used
length(meth)

meth["Sex"] <- ""
meth[grep('ACE_score_extended', names(meth), value = T)] <-
  paste0('~I(', paste0('as.integer(', ACEs, ')', collapse = '+'), ')')
meth[grep('ACE_score_classic', names(meth), value = T)] <-
  paste0('~I(',
         paste0('as.integer(', classicACEs, ')', collapse = '+'),
         ')')
meth["UI"] <-
  '~ I(ifelse(Daywetting == "No" & Bedwetting =="No", "No", "Yes"))' # passive impute

## Checks lengths [should be 100 variables]
length(meth)
dim(pred_mat)
table(rowSums(pred_mat)) # nr predictors per variable
names(which(rowSums(pred_mat) == 0)) # complete variables

# Save the datasets
save(
  meth,
  pred_mat,
  dat_boys,
  dat_girls,
  file = paste0(loc_out, '04-init-imputation-data.RData')
)

#===============================================================================
# Code block 8: Look at missing data again for some reason
#===============================================================================
miss_boys <-
  unlist(lapply(dat_boys, function(x)
    sum(is.na(x)))) / nrow(dat_boys)
sort(miss_boys[miss_boys > 0.25], decreasing = TRUE)

miss_girls <-
  unlist(lapply(dat_girls, function(x)
    sum(is.na(x)))) / nrow(dat_girls)
sort(miss_girls[miss_girls > 0.25], decreasing = TRUE)


#===============================================================================
# Code block 9: Run Imputations
# Run for males and females
#===============================================================================
catf('\n\n',
     paste0(Sys.time()),
     'Code block 4: Imputation-RECOMMEND HPC\n')

#### ### Run imputations tests
#### start <- Sys.time()
#### imp_boys <- mice(
####   dat_boys,
####   m = 4,
####   maxit = 1,
####   print = TRUE,
####   method = meth,
####   predictorMatrix = pred_mat,
####   seed = 140817
#### )
#### print(Sys.time() - start)
#### # Seq time = 1.19 mins
#### 
#### # Parlmice
#### start <- Sys.time()
#### imp_boys2 <- parlmice(
####   dat_boys,
####   m = 4,
####   maxit = 1,
####   print = TRUE,
####   method = meth,
####   predictorMatrix = pred_mat,
####   cluster.seed=140817,
####   n.core = 4,
####   n.imp.core = 1,
####   cl.type = "FORK"
#### )
#### print(Sys.time() - start)
#### # Seq time = 31.51 seconds
#### 
#### # Parlmice
#### start <- Sys.time()
#### imp_boys3 <- parlmice(
####   dat_boys,
####   m = 4,
####   maxit = 1,
####   print = TRUE,
####   method = meth,
####   predictorMatrix = pred_mat,
####   cluster.seed=140817,
####   n.core = 4,
####   n.imp.core = 1,
####   cl.type = "FORK"
#### )
#### print(Sys.time() - start)
#### # Seq time = 31.51 seconds
#### 
#### # the above work with logreg_2 and parallel is indeed parallel
#### 
#### summary(imp_boys[["imp"]][["Physical_abuse"]])
#### summary(imp_boys2[["imp"]][["Physical_abuse"]])
#### summary(imp_boys3[["imp"]][["Physical_abuse"]])
#### # Expect that 1 and 2/3 will differ = TRUE
#### # Expect that 2 and 3 be the same = TRUE

# Parlmice
start <- Sys.time()
imp_boys <- parlmice(
  dat_boys,
  m = 80,
  maxit = 50,
  print = TRUE,
  method = meth,
  predictorMatrix = pred_mat,
  cluster.seed=140817,
  n.core = 20,
  n.imp.core = 4,
  cl.type = "FORK"
)
print(Sys.time() - start)

start <- Sys.time()
imp_girls <- parlmice(
  dat_girls,
  m = 80,
  maxit = 50,
  print = TRUE,
  method = meth,
  predictorMatrix = pred_mat,
  cluster.seed=140817,
  n.core = 20,
  n.imp.core = 4,
  cl.type = "FORK"
)
print(Sys.time() - start)

# Save the datasets
save(
  imp_boys,
  imp_girls,
  file = paste0(loc_out, '05-imputation-data.RData')
)

#===============================================================================
# Code block 10: Plot diagnostics
# Stripplots, convergence checks, box plots, proportions, density plots
#===============================================================================

catf('\n\n',paste0(Sys.time()),
     'Code block 5: Check imputation\n')

# Change logicals to factors
imp_boys_sub <-
  complete(imp_boys, 'long', include = TRUE) %>% mutate(across(
    c(
      "Physical_abuse",
      "Sexual_abuse",
      "Emotional_abuse",
      "Emotional_neglect",
      "Bullying",
      "Violence_parents",
      "Substance_abuse",
      "Mental_hlth_suicide",
      "Convicted_offences",
      "Parental_separation",
      "Social_class",
      "Financial_difficulties",
      "Neighbourhood",
      "Social_support_parents",
      "Physical_illness_child",
      "Physical_illness_parents",
      "Parent_child_bond",
      "Edu_partner_mthrReport",
      "Edu_ptnr_self",
      "Edu_mother_ptnrReport",
      "Ptnr_emo_cruel_18yrs",
      "Antidep_use_mthr_18yrs",
      "PtnrChild_phys_force_18yrs",
      "PtnrChild_sev_phys_force_18yrs",
      "PtnrChild_kiss_touch_18yrs",
      "PtnrChild_sev_kiss_touch_18yrs",
      "PtnrChild_pressure_sex_18yrs",
      "PtnrChild_physforce_sex_18yrs",
      "PtnrChild_feel_scared_18yrs"
    ),
    as_factor
  ))
imp_boys_sub <- as.mids(imp_boys_sub)

imp_girls_sub <-
  complete(imp_girls, 'long', include = TRUE) %>% mutate(across(
    c(
      "Physical_abuse",
      "Sexual_abuse",
      "Emotional_abuse",
      "Emotional_neglect",
      "Bullying",
      "Violence_parents",
      "Substance_abuse",
      "Mental_hlth_suicide",
      "Convicted_offences",
      "Parental_separation",
      "Social_class",
      "Financial_difficulties",
      "Neighbourhood",
      "Social_support_parents",
      "Physical_illness_child",
      "Physical_illness_parents",
      "Parent_child_bond",
      "Edu_partner_mthrReport",
      "Edu_ptnr_self",
      "Edu_mother_ptnrReport",
      "Ptnr_emo_cruel_18yrs",
      "Antidep_use_mthr_18yrs",
      "PtnrChild_phys_force_18yrs",
      "PtnrChild_sev_phys_force_18yrs",
      "PtnrChild_kiss_touch_18yrs",
      "PtnrChild_sev_kiss_touch_18yrs",
      "PtnrChild_pressure_sex_18yrs",
      "PtnrChild_physforce_sex_18yrs",
      "PtnrChild_feel_scared_18yrs"
    ),
    as_factor
  ))
imp_girls_sub <- as.mids(imp_girls_sub)

### Combine genders
imp_all<-rbind(imp_boys,imp_girls)
save(imp_all,file='06-imp-all.RData')

imp_all_sub<-rbind(imp_boys_sub,imp_girls_sub)
save(imp_all_sub,file='06b-imp-all-sub.RData')

### Check convergence of iterations
pdf('02-check-conv-boys.pdf')
plot(imp_boys)
dev.off()

pdf('03-check-conv-girls.pdf')
plot(imp_girls)
dev.off()

### Load merged imp data
#load("06-imp_all.RData")

# Box plot of ACE score
pdf('04-box_plot_ace_score.pdf')
bwplot(imp_all_sub, ACE_score_classic ~.imp)
dev.off()

# Density and box plots of continuous variables
pdf('05-density_plot.pdf')
densityplot(imp_all_sub, ~ACE_score_classic + BMI + SDQ + Parity + Birthweight, par.strip.text = list(cex=0.75, font=1), main="Density of observed (blue) and imputed (red) values")
dev.off()

pdf('06-box_plot_continuous.pdf')
bwplot(imp_all_sub, ACE_score_classic + BMI + SDQ + Parity + Birthweight ~.imp)
dev.off()

# Create binary and categorical datasets for proportional plots
source("../scripts/02-ACEs/03e_prop_plot_bin.R")
source("../scripts/02-ACEs/03f_prop_plot_cat.R")

imp_all_sub_binary <- complete(imp_all_sub, 'long', include = TRUE) %>% select(!c(Edu_ptnr_self, Edu_mother_ptnrReport,
                                                                        Edu_partner_mthrReport, mateducation, House_tenure,
                                                                        Marital_status, Crowding, smoked_pregnancy_cat))
imp_all_sub_binary <- as.mids(imp_all_sub_binary)


imp_all_sub_cat <- complete(imp_all_sub, 'long', include = TRUE) %>% select(.imp, .id, Edu_ptnr_self, Edu_mother_ptnrReport,
                                                                         Edu_partner_mthrReport, mateducation, House_tenure,
                                                                         Marital_status, Crowding, smoked_pregnancy_cat)
imp_all_sub_cat <- as.mids(imp_all_sub_cat)

### Proportional plots
pdf("07-proportions_binary.pdf")
propplot_bin(imp_all_sub_binary)
dev.off()

pdf("08-proportions_cats.pdf")
propplot_cat(imp_all_sub_cat)
dev.off()

imp_all_sub_binary_b <- complete(imp_all_sub, 'long', include = TRUE) %>% select(.imp, .id, UI, Daywetting, Bedwetting, Nocturia, Urgency, Frequency, Voiding_postponement, Voiding_volume)
imp_all_sub_binary_b <- as.mids(imp_all_sub_binary_b)

### Proportional plots
pdf("07-proportions_binary_b.pdf")
propplot_bin(imp_all_sub_binary_b)
dev.off()

#===============================================================================
# Code block 11: Checks and formatting
# From Lottes code
#===============================================================================
# Change logicals to factors
imp_all_sub <-
  complete(imp_all, 'long', include = TRUE) %>% mutate(across(
    c(
      "Physical_abuse",
      "Sexual_abuse",
      "Emotional_abuse",
      "Emotional_neglect",
      "Bullying",
      "Violence_parents",
      "Substance_abuse",
      "Mental_hlth_suicide",
      "Convicted_offences",
      "Parental_separation",
      "Social_class",
      "Financial_difficulties",
      "Neighbourhood",
      "Social_support_parents",
      "Physical_illness_child",
      "Physical_illness_parents",
      "Parent_child_bond",
      "Edu_partner_mthrReport",
      "Edu_ptnr_self",
      "Edu_mother_ptnrReport",
      "Ptnr_emo_cruel_18yrs",
      "Antidep_use_mthr_18yrs",
      "PtnrChild_phys_force_18yrs",
      "PtnrChild_sev_phys_force_18yrs",
      "PtnrChild_kiss_touch_18yrs",
      "PtnrChild_sev_kiss_touch_18yrs",
      "PtnrChild_pressure_sex_18yrs",
      "PtnrChild_physforce_sex_18yrs",
      "PtnrChild_feel_scared_18yrs"
    ),
    as_factor
  ))
imp_all_sub <- as.mids(imp_all_sub)

imp_all_sub2 <-
  complete(imp_all_sub, 'long', include = TRUE) %>% mutate(across(
    c(
      "Physical_abuse",
      "Sexual_abuse",
      "Emotional_abuse",
      "Emotional_neglect",
      "Bullying",
      "Violence_parents",
      "Substance_abuse",
      "Mental_hlth_suicide",
      "Convicted_offences",
      "Parental_separation"
    ),
    ~ case_match(., "0" ~ 0, "1" ~ 1)
  ))
imp_all_sub2 <- as.mids(imp_all_sub2)

com_all <- complete(imp_all_sub2, "long", include = TRUE)
table(is.na(com_all$ACEcat_extended_0_8))
nrow(com_all[com_all$.imp==0,]);colSums(is.na(com_all[com_all$.imp==0,]))
write.csv(com_all,file=paste0(loc_inp,'07-com-all.csv'))

description_imputation<-
  data.frame('original_var'=names(imp_all_sub$method),
             'Variable'=gsub('_org|_dup','',names(imp_all_sub$method)),
             'Type of variable'=sapply(
               imp_all_sub$data[,match(names(imp_all_sub$method),colnames(imp_all_sub$data))],
               function(x) {cl=class(x);paste0(cl,
                                               ifelse(cl%in%c("character", "factor"),
                                                      paste0(' (',length(levels(factor(x))),')'),''))}),
             'Regression model to predict missing in this variable'=imp_all_sub$method)


description_imputation[[3]]<-gsub("factor (2)",'dichotomous', fixed = TRUE ,description_imputation[[3]])
description_imputation[[3]]<-gsub('integer','continuous',description_imputation[[3]])
description_imputation[[3]]<-gsub('numeric','continuous',description_imputation[[3]])
description_imputation[[3]]<-gsub('factor','categorical',description_imputation[[3]])
description_imputation[[4]]<-gsub('logreg','Logistic regression',description_imputation[[4]])
description_imputation[[4]]<-gsub('pmm','Predictive mean matching',description_imputation[[4]])
description_imputation[[4]]<-gsub('polyreg','Polytomous (unordered) regression',description_imputation[[4]])

description_imputation$Variable<-gsub('_',' ',description_imputation$Variable)
description_imputation$Variable<-gsub('0 8yrs','0-8yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Parent child bond','Parent-child bond',description_imputation$Variable)
description_imputation$Variable<-gsub('ACE score classic','ACE-score',description_imputation$Variable)
description_imputation$Variable<-gsub('Violence parents','Violence between parents',description_imputation$Variable)
description_imputation$Variable<-gsub('Substance abuse','Household substance abuse',description_imputation$Variable)
description_imputation$Variable<-gsub('Mental hlth suicide','Mental health problems or suicide',description_imputation$Variable)
description_imputation$Variable<-gsub('Convicted offences','Parents convicted of criminal offence',description_imputation$Variable)
description_imputation$Variable<-gsub('Parental separation','Parental separation',description_imputation$Variable)

description_imputation$Variable<-gsub('sc household 18wgest','Household social class at 18wks gestation',description_imputation$Variable)
description_imputation$Variable<-gsub('Edu partner mthrReport','Mother-reported highest educational level partner',description_imputation$Variable)
description_imputation$Variable<-gsub('Edu ptnr self','Self-reported highest educational level partner ',description_imputation$Variable)
description_imputation$Variable<-gsub('Edu mother ptnrReport','Partner-reported highest educational level mother',description_imputation$Variable)
description_imputation$Variable<-gsub('EPDS mthr 18wksG','Maternal depression score (EPDS) at 18 wks gestation',description_imputation$Variable)
description_imputation$Variable<-gsub('EPDS mthr 32wksG','Maternal depression score (EPDS) at 32 wks gestation',description_imputation$Variable)
description_imputation$Variable<-gsub('EPDS ptnr 18wksG','Partner depression score (EPDS) at 18 wks gestation',description_imputation$Variable)
description_imputation$Variable<-gsub('Maternal weight','Maternal pre-pregnancy weight (Kg)',description_imputation$Variable)
description_imputation$Variable<-gsub('Maternal BMI','Maternal pre-pregnancy BMI',description_imputation$Variable)
description_imputation$Variable<-gsub('AUDIT ptnr 18yrs','AUDIT Alcohol disorders score (partner)',description_imputation$Variable)
description_imputation$Variable<-gsub('AUDIT mthr 18yrs','AUDIT Alcohol disorders score (mother)',description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild phys force 18yrs',"YPs age when partners have used physical force",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild kiss touch 18yrs',"YPs age when YP's partners have pressured them into kissing/touching/something else",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild physforce sex 18yrs',"YPs age when partners have pressured them into having sexual intercourse (severe)",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild sev phys force 18yrs',"YPs age when partners have used physical force (severe)",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild sev kiss touch 18yrs',"YPs age when YP's partners have pressured them into kissing/touching/something else (severe)",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild pressure sex 18yrs',"YPs age when partners have pressured them into having sexual intercourse",description_imputation$Variable)
description_imputation$Variable<-gsub('PtnrChild feel scared 18yrs',"YPs age when partners behaviour  made them feel scared or frightened",description_imputation$Variable)
description_imputation$Variable<-gsub('Ptnr emo cruel 18yrs',"Respondent's partner was emotionally cruel to respondent's children",description_imputation$Variable)
description_imputation$Variable<-gsub('Antidep use mthr 18yrs',"Mother's antidepressant use",description_imputation$Variable)

description_imputation$Variable<-gsub('Sex','Sex',description_imputation$Variable)
description_imputation$Variable<-gsub('Ethnicity','Ethnicity of child',description_imputation$Variable)
description_imputation$Variable<-gsub('Birthweight','Birthweight of child in grams',description_imputation$Variable)
description_imputation$Variable<-gsub('Gestational age','Gestational age in weeks at delivery',description_imputation$Variable)
description_imputation$Variable<-gsub('Parity','Parity',description_imputation$Variable)
description_imputation$Variable<-gsub('matage','Maternal age in years at delivery',description_imputation$Variable)
description_imputation$Variable<-gsub('mateducation','Maternal education (pregnancy)',description_imputation$Variable)
description_imputation$Variable<-gsub('House tenure','Housing tenure (antenatal)',description_imputation$Variable)
description_imputation$Variable<-gsub('devlevel','Childs developmental level (score) at age 18 months ',description_imputation$Variable)
description_imputation$Variable<-gsub('matage','Maternal age in years at delivery',description_imputation$Variable)
description_imputation$Variable<-gsub('Marital status','Marital status of mother (antenatal)',description_imputation$Variable)
description_imputation$Variable<-gsub('Crowding','Home crowding index (antenatal)',description_imputation$Variable)
description_imputation$Variable<-gsub('smoked pregnancy cat','Maternal smoking during pregnancy',description_imputation$Variable)
description_imputation$Variable<-gsub('age bloods 9','Age of child at blood draw',description_imputation$Variable)

description_imputation$Variable<-gsub('Daywetting','Daytime wetting at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Bedwetting','Bedtime wetting at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Frequency','Frequent urination at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Nocturia','Nocturia at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Urgency','Urgency for urination at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Voiding postponement','Voiding postponement at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('Voiding volume','High voiding volume at age 14 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('UI','Any urinary incontinence at age 14 yrs',description_imputation$Variable)

description_imputation$Variable<-gsub('daywet 7 9','Daytime wetting from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('bedwet 7 9','Bedtime wetting from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('freq urine 7 9','Frequent urination from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('nocturia 7 9','Nocturia from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('urgency 7 9','Urgency for urination from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('void postpone 7 9','Voiding postponement from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('signs 7 9','Signs of urination from 7-9 yrs',description_imputation$Variable)
description_imputation$Variable<-gsub('day soiling 7 9','Daytime soiling from 7-9 yrs',description_imputation$Variable)

write.csv(description_imputation,file=paste0(loc_inp,'08-description-imputation.csv'))

##Calculate spread in imputed datasets (cc vs imp mean, imp sd, imp min, imp max):
calc_ACE_prev<-function(i_dataset,imp_long,ACEs,...){
  require(matrixStats)
  colSums(sapply(imp_long[imp_long$.imp%in%i_dataset,ACEs],as.numeric),na.rm=T)/
    colSums(!is.na(imp_long[imp_long$.imp%in%i_dataset,ACEs]))
}

ACEs <- c(
  "Physical_abuse",
  "Sexual_abuse",
  "Emotional_abuse",
  "Emotional_neglect",
  "Bullying",
  "Violence_parents",
  "Substance_abuse",
  "Mental_hlth_suicide",
  "Convicted_offences",
  "Parental_separation"
)

ACE_prevalence_per_dataset<-sapply(as.character(unique(com_all$.imp)),calc_ACE_prev,imp_long=com_all,ACEs=ACEs)

SItab5_compare_imputation_all<-
  data.frame(#'Nr_imputed_datasets'=rep(ncol(ACE_prevalence_per_dataset)-1,nrow(ACE_prevalence_per_dataset,)),
    'ACE'=ACEs,
    'Complete cases (n)'=colSums(!is.na(com_all[com_all$.imp%in%0,ACEs])),
    'Missingness (%)'=round(colSums(is.na(com_all[com_all$.imp%in%0,ACEs]))/
                              nrow(com_all[com_all$.imp%in%0,ACEs])*100,2),
    'Prevalence complete cases (%)'=round(ACE_prevalence_per_dataset[,1]*100,2),
    'Prevalence imputed data (mean %)'=
      round(rowMeans(ACE_prevalence_per_dataset[,2:ncol(ACE_prevalence_per_dataset),drop=F])*100,2),
    'Fold increase in prevalence from complete cases to mean prevalence in imputed data'=
      round(
        (rowMeans(ACE_prevalence_per_dataset[,2:ncol(ACE_prevalence_per_dataset),drop=F])-ACE_prevalence_per_dataset[,1])/
          ACE_prevalence_per_dataset[,1],2),
    'Lowest prevalence across all imputed datasets'=
      round(rowMins(ACE_prevalence_per_dataset[,2:ncol(ACE_prevalence_per_dataset),drop=F])*100,2),
    'Highest prevalence across all imputed datasets'=
      round(rowMaxs(ACE_prevalence_per_dataset[,2:ncol(ACE_prevalence_per_dataset),drop=F])*100,2),
    check.names=F)

rownames(SItab5_compare_imputation_all)<-gsub('_',' ',rownames(SItab5_compare_imputation_all))
rownames(SItab5_compare_imputation_all)<-gsub('0 8yrs','0-8yrs',rownames(SItab5_compare_imputation_all))
rownames(SItab5_compare_imputation_all)<-gsub('Parent child bond','parent-child bond',
                                              rownames(SItab5_compare_imputation_all))

write.csv(SItab5_compare_imputation_all,file=paste0(loc_out,'09-compare-imputation-all.csv'),row.names = F)

################
#End of script
################
