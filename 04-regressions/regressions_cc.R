# Mediation of inflammation in the associations between ACEs and LUTS
#===============================================================================
# Kimberley Burrows
# Regressions for Table S9, S10 and S8b(omitted) using complete case data
#===============================================================================

#-------------------------------------------------------------------------------
# Admin ------------------------------------------------------------------------
rm(list = ls()) #clear global R environment
renv::status()
renv::snapshot()

# Library
library(mice)
library(kableExtra)
library(modelsummary)
library(miceadds)
library(flextable)
library(tidyverse)
library(gtsummary)
library(oddsratio)
library(gt)
library(forcats)


# Notes: -----------------------------------------------------------------------
# Paper Table S9: LUTS ~ ACE score [+ confs]
# Paper Table S10: LUTS ~ Inflammation [+ conf + ACE score + L's]
# Paper Table S8b (omitted): Inflammation ~ ACEs [+ confs] & Inflammation ~ Inflammation [+ confs + ACE score + L's]
# Date in: imputation dataset from mice(): class mids
#-------------------------------------------------------------------------------

#===============================================================================
# Options and data prep
#===============================================================================

## Define paths for input/output
Tables <-
  "B:/Kim/02-Inflammation-mediation/02-analysis/tables/"
Figures <-
  "B:/Kim/02-Inflammation-mediation/02-analysis/figures/"
Data <-
  "B:/Kim/02-Inflammation-mediation/02-analysis/data/"
Scripts <-
  "B:/Kim/02-Inflammation-mediation/02-analysis/scripts/"

setwd(Data)

## Load data and convert to long
load("HPC/06b-imp-all-sub.RData")
imps_long <- complete(imp_all_sub, action = 'long', include = TRUE)

## rename some vars
names(imps_long)[names(imps_long) == "ACE_score_classic"] <- "ACE_Score"
names(imps_long)[names(imps_long) == "devlevel"] <- "Developmental_stage"
names(imps_long)[names(imps_long) == "matage"] <- "Maternal_age"
names(imps_long)[names(imps_long) == "smoked_pregnancy_cat"] <- "Smoking"
names(imps_long)[names(imps_long) == "mateducation"] <- "Maternal_education"

## Select model vars
imps_long <- imps_long %>% select(
  .imp,
  .id,
  Physical_abuse,
  Sexual_abuse,
  Emotional_abuse,
  Emotional_neglect,
  Bullying,
  Violence_parents,
  Substance_abuse,
  Mental_hlth_suicide,
  Convicted_offences,
  Parental_separation,
  ACE_Score,
  IL6,
  CRP,
  UI,
  Daywetting,
  Bedwetting,
  Nocturia,
  Frequency,
  Urgency,
  Voiding_postponement,
  Voiding_volume,
  BMI,
  SDQ,
  Constipation,
  Sex,
  Birthweight,
  Gestational_age,
  Ethnicity,
  Developmental_stage,
  Maternal_age,
  Smoking,
  Parity,
  Maternal_education,
  Marital_status,
  House_tenure,
  Crowding,
  Age_blood_draw_9
) # 39 vars

## Generate log variables
imps_long$logIL6 <- log(imps_long$IL6)
imps_long$logCRP <- log(imps_long$CRP)

imps_long$logIL6 <- imps_long$IL6
imps_long$logCRP <- imps_long$CRP

## Collapse ACE score 7+ to mitigate sparse cells
imps_long$ACE_Score_cat <-as.factor(imps_long$ACE_Score)
imps_long <- imps_long %>%
  mutate(ACE_score = fct_collapse(ACE_Score_cat, "7"=c("7", "8", "9", "10")))
imps_long$ACE_score <- as.numeric(as.character(imps_long$ACE_score))

## recode the factors (FALSE/TRUE to No/Yes)
imps_long[, c(3:12)] <- lapply(imps_long[, c(3:12)], function(x) plyr::revalue(x, c("0"="No", "1"="Yes")))
imp_use <- as.mids(imps_long)

## save version of this imputation dataset for future use
#saveRDS(imp_use, file = "imps_to_use.RData")

## Keep the complete case data
complete_cases <- imp_use$data %>%
  filter(complete.cases(.))
#===============================================================================

#===============================================================================
## TABLE 1. LUTS ~ ACEs [+ conf]

conf <- c("Sex",
          "Birthweight",
          "Gestational_age",
          "Ethnicity",
          "Developmental_stage",
          "Maternal_age",
          "Smoking",
          "Parity",
          "Maternal_education",
          "Marital_status",
          "House_tenure",
          "Crowding")

# vectors of model variables
vars <-
  c("ACE_score", paste("ACE_score + ", paste(conf, collapse = " + ")))
outcomes <- c("UI", "Daywetting", "Bedwetting", "Urgency", "Nocturia", "Frequency", "Voiding_postponement", "Voiding_volume")

# tables of results for each model unadjusted vs. adjusted (8 outcomes)
mods <- list()
for (j in seq_along(outcomes)) {
  mods[[j]] <- list()
  for (i in seq_along(vars)) {
    mods[[j]][[i]] <- glm(as.formula(paste0(outcomes[[j]], " ~ ", vars[[i]])), family = binomial(link = "logit"), data = complete_cases) %>%
      tbl_regression(
        exponentiate = TRUE,
        pvalue_fun = ~style_sigfig(., digits = 3),
        show_single_row = ACE_score,
        label = list(ACE_score ~ outcomes[[j]])
      ) %>%
      #bold_p() %>%
      modify_table_body(filter,
                        !(variable %in% conf)) %>%
      modify_header(label = "**ACE score on LUTS**") %>%
      modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
  }}

results <- lapply(1:8, function(k) {
  tbl_merge(mods[[k]],
            tab_spanner = c("**Univariable**", "**Multivariable**")) %>%
    modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
})

table <- tbl_stack(list(results[[1]],
               results[[2]],
               results[[3]],
               results[[4]],
               results[[5]],
               results[[6]],
               results[[7]],
               results[[8]]
))

table <- table %>%
  as_gt()

gtsave(table, file = paste0(Tables,"table_1_cc_s9.docx"))

rm(table, results, mods)
gc()
#===============================================================================

#===============================================================================
## TABLE 2. LUTS ~ Inflammation [+ conf + ACEs + L's]
conf2 <- c("Sex",
          "Birthweight",
          "Gestational_age",
          "Ethnicity",
          "Developmental_stage",
          "Maternal_age",
          "Smoking",
          "Parity",
          "Maternal_education",
          "Marital_status",
          "House_tenure",
          "Crowding",
          "BMI",
          "SDQ",
          "Age_blood_draw_9",
          "Constipation",
          "ACE_score")

# vectors of model variables
## IL6 ##
vars <-
  c("logIL6", paste("logIL6 + ", paste(conf2, collapse = " + ")))
outcomes <- c("UI", "Daywetting", "Bedwetting", "Urgency", "Nocturia", "Frequency", "Voiding_postponement", "Voiding_volume")

# tables of results for each model unadjusted vs. adjusted
mods <- list()
for (j in seq_along(outcomes)) {
  mods[[j]] <- list()
  for (i in seq_along(vars)) {
    mods[[j]][[i]] <- glm(as.formula(paste0(outcomes[[j]], " ~ ", vars[[i]])), family = binomial(link = "logit"), data = complete_cases) %>%
      tbl_regression(
        exponentiate = TRUE,
        pvalue_fun = ~style_sigfig(., digits = 3),
        show_single_row = logIL6,
        label = list(logIL6 ~ outcomes[[j]])
      ) %>%
      modify_table_body(filter,
                        !(variable %in% conf2)) %>%
      modify_header(label = "**Outcome**") %>%
      modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
  }}

results_il6 <- lapply(1:8, function(k) {
  tbl_merge(mods[[k]],
            tab_spanner = c("**Univariable**", "**Multivariable**")) %>%
    modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
})

## CRP ##
vars <-
  c("logCRP", paste("logCRP + ", paste(conf2, collapse = " + ")))
outcomes <- c("UI", "Daywetting", "Bedwetting", "Urgency", "Nocturia", "Frequency", "Voiding_postponement", "Voiding_volume")

# tables of results for each model unadjusted vs. adjusted
mods <- list()
for (j in seq_along(outcomes)) {
  mods[[j]] <- list()
  for (i in seq_along(vars)) {
    mods[[j]][[i]] <- glm(as.formula(paste0(outcomes[[j]], " ~ ", vars[[i]])), family = binomial(link = "logit"), data = complete_cases) %>%
      tbl_regression(
        exponentiate = TRUE,
        pvalue_fun = ~style_sigfig(., digits = 3),
        show_single_row = logCRP,
        label = list(logCRP ~ outcomes[[j]])
      ) %>%
      modify_table_body(filter,
                        !(variable %in% conf2)) %>%
      modify_header(label = "**Outcome**") %>%
      modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
  }}

results_crp <- lapply(1:8, function(k) {
  tbl_merge(mods[[k]],
            tab_spanner = c("**Univariable**", "**Multivariable**")) %>%
    modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
})

# stack IL6 and CRP results
table <- tbl_stack(list(results_il6[[1]],
                          results_il6[[2]],
                          results_il6[[3]],
                          results_il6[[4]],
                          results_il6[[5]],
                          results_il6[[6]],
                          results_il6[[7]],
                          results_il6[[8]],
                          results_crp[[1]],
                          results_crp[[2]],
                          results_crp[[3]],
                          results_crp[[4]],
                          results_crp[[5]],
                          results_crp[[6]],
                          results_crp[[7]],
                          results_crp[[8]]
)) %>%
  as_gt() %>%
  #fmt_number(decimals =2) %>%
  tab_row_group(
  label = "logIL6 (pg/ml)",
  rows = 1:8
) %>%
  tab_row_group(
    label = "logCRP (mg/L)",
    rows = 9:16
  ) %>%
  row_group_order(groups = c("logIL6 (pg/ml)", "logCRP (mg/L)")) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgrey"),
      cell_text(style = "italic")
    ),
    locations = cells_row_groups()
  ) %>%
  tab_options(row_group.padding = px(0.8)) %>%
  cols_align(
    align = "left",
    columns = everything()
  )
table

gtsave(table, file = paste0(Tables,"table_2_cc_s10.docx"))

# P values - need to find a way to format properly in above code, for now add in manually
results_il6[[1]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[2]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[3]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[4]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[5]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[6]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[7]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[8]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_il6[[1]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[2]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[3]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[4]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[5]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[6]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[7]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_il6[[8]][["tbls"]][[2]][["table_body"]][["p.value"]]

results_crp[[1]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[2]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[3]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[4]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[5]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[6]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[7]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[8]][["tbls"]][[1]][["table_body"]][["p.value"]]
results_crp[[1]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[2]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[3]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[4]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[5]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[6]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[7]][["tbls"]][[2]][["table_body"]][["p.value"]]
results_crp[[8]][["tbls"]][[2]][["table_body"]][["p.value"]]

rm(table, results_crp, results_il6, mods)
gc()
#===============================================================================

#===============================================================================
## TABLE S8b. Inflammation ~ ACEs [+ confs] - Note that this was not included in the supplementary for some reason
conf <- c("Sex",
          "Birthweight",
          "Gestational_age",
          "Ethnicity",
          "Developmental_stage",
          "Maternal_age",
          "Smoking",
          "Parity",
          "Maternal_education",
          "Marital_status",
          "House_tenure",
          "Crowding")


# vectors of model variables
vars <-
  c("ACE_score", paste("ACE_score + ", paste(conf, collapse = " + ")))
outcomes <- c("logIL6", "logCRP")

# tables of results for each model unadjusted vs. adjusted
mods <- list()
for (j in seq_along(outcomes)) {
  mods[[j]] <- list()
  for (i in seq_along(vars)) {
    mods[[j]][[i]] <- glm(as.formula(paste0(outcomes[[j]], " ~ ", vars[[i]])), family = gaussian(), data = complete_cases) %>%
      tbl_regression(
        exponentiate = TRUE,
        pvalue_fun = ~style_sigfig(., digits = 3),
        show_single_row = ACE_score,
        label = list(ACE_score ~ outcomes[[j]])
      ) %>%
      #bold_p() %>%
      modify_table_body(filter,
                        !(variable %in% conf)) %>%
      modify_header(label = "**Outcome**") %>%
      modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
  }}

results <- lapply(1:2, function(k) {
  tbl_merge(mods[[k]],
            tab_spanner = c("**Univariable**", "**Multivariable**")) %>%
    modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
})

# Now CRP ~ IL6 + conf + ACE score + L's
conf2 <- c("Sex",
           "Birthweight",
           "Gestational_age",
           "Ethnicity",
           "Developmental_stage",
           "Maternal_age",
           "Smoking",
           "Parity",
           "Maternal_education",
           "Marital_status",
           "House_tenure",
           "Crowding",
           "BMI",
           "SDQ",
           "Age_blood_draw_9",
           "Constipation",
           "ACE_score")

vars <-
  c("logIL6", paste("logIL6 + ", paste(conf2, collapse = " + ")))
outcomes <- c("logCRP")

# tables of results for each model unadjusted vs. adjusted
mods <- list()
for (j in seq_along(outcomes)) {
  mods[[j]] <- list()
  for (i in seq_along(vars)) {
    mods[[j]][[i]] <- glm(as.formula(paste0(outcomes[[j]], " ~ ", vars[[i]])), family = gaussian(), data = complete_cases) %>%
      tbl_regression(
        pvalue_fun = ~style_sigfig(., digits = 3),
        exponentiate = TRUE#,
        #show_single_row = ACE_score,
        #label = list(ACE_score ~ outcomes[[j]])
      ) %>%
      #bold_p() %>%
      modify_table_body(filter,
                        !(variable %in% conf2)) %>%
      modify_header(label = "**Outcome**") %>%
      modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
  }}

results2 <- lapply(1:1, function(k) {
  tbl_merge(mods[[k]],
            tab_spanner = c("**Univariable**", "**Multivariable**")) %>%
    modify_footnote(update = everything() ~ NA, abbreviation = TRUE)
})

# Stack tables
table <- tbl_stack(list(results[[1]],
                        results[[2]],
                        results2[[1]]
)) %>%
  as_gt() %>%
  #fmt_number(decimals =2) %>%
  tab_row_group(
    label = "ACE Score",
    rows = 1:2
  ) %>%
  tab_row_group(
    label = "logIL6 (pg/ml)",
    rows = 3
  ) %>%
  row_group_order(groups = c("ACE Score", "logIL6 (pg/ml)")) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgrey"),
      cell_text(style = "italic")
    ),
    locations = cells_row_groups()
  ) %>%
  tab_options(row_group.padding = px(0.5)) %>%
  cols_align(
    align = "left",
    columns = everything()
  )

table

gtsave(table, file = paste0(Tables,"table_s8_cc.docx"))

rm(table, results, results2, mods, conf2)
gc()

# END