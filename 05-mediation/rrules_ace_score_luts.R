# Mediation of inflammation in the associations between ACEs and LUTS
#===============================================================================
# Kimberley Burrows
# Completed:
#===============================================================================

# Notes:
# Code adapted from Andrea Cortes: 6_ACEs_PEs_extracting info from log.R
# Additions: file appending, Rubins Rules, and tables
# input files are Stata log files:
## g-formula logs with "return" called after model has run to give list of scalars with
# model estimates needed.
# HPC used to run each imputation set individually
# Log files therefore need appending


#-------------------------------------------------------------------------------
# Admin ------------------------------------------------------------------------
rm(list = ls())#clear global R environment
#renv::status()
#renv::snapshot()

# Library
library(data.table) #v.1.14.0
library(readr)
library(stringr)
library(reshape2)
library(tidyverse)
library(gt)

setwd("B:/Kim/02-Inflammation-mediation/02-analysis/scripts/05-mediation")

#===============================================================================
# Options and data prep
# Append the individual mediation analysis logs
# Easier done in shell:
## cat *.txt >> all_logs.txt
# Note that the order of files appended follows convention of
# 1, 11, 12 ... 2, 21, 22... due to lack of leading zeros in file name
# Does not matter for RRs
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

log_ui <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_UI.txt")
log_bw <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Bedwetting.txt")
log_dw <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Daywetting.txt")
log_freq <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Frequency.txt")
log_noct <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Nocturia.txt")
log_urgency <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Urgency.txt")
log_void_post <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Voiding_postponement.txt")
log_void_vol <- readLines("B:/Kim/02-Inflammation-mediation/02-analysis\\data\\HPC\\mediation_logs\\all_logs_ACE_score_Voiding_volume.txt")

#===============================================================================
# Extract estimates from log files:
#===============================================================================

# any UI
DT_log <- data.table(txt = log_ui[grepl(pattern = 'Imp number is:  |) =', log_ui)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(ui_data <- gformula_dataset[, col_order])

# Daytime wetting
DT_log <- data.table(txt = log_dw[grepl(pattern = 'Imp number is:  |) =', log_dw)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(dw_data <- gformula_dataset[, col_order])

# Bedtime wetting
DT_log <- data.table(txt = log_bw[grepl(pattern = 'Imp number is:  |) =', log_bw)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(bw_data <- gformula_dataset[, col_order])

# Frequency
DT_log <- data.table(txt = log_freq[grepl(pattern = 'Imp number is:  |) =', log_freq)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(freq_data <- gformula_dataset[, col_order])

# Nocturia
DT_log <- data.table(txt = log_noct[grepl(pattern = 'Imp number is:  |) =', log_noct)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(noct_data <- gformula_dataset[, col_order])

# Urgency
DT_log <- data.table(txt = log_urgency[grepl(pattern = 'Imp number is:  |) =', log_urgency)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(urgency_data <- gformula_dataset[, col_order])

# Voiding postponement
DT_log <- data.table(txt = log_void_post[grepl(pattern = 'Imp number is:  |) =', log_void_post)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(void_post_data <- gformula_dataset[, col_order])

# Voiding volume
DT_log <- data.table(txt = log_void_vol[grepl(pattern = 'Imp number is:  |) =', log_void_vol)]) #to convert it to 1 column data table + get only those columns with the info I need.
DT_log <- DT_log[!grepl(". di ", DT_log$txt),]
DT_log[, grp := cumsum(grepl('Imp number is', txt))] #putting a number per each imputation, based on when MC_sims appears.
str(DT_log)
DT_log$var2 <- parse_number(DT_log$txt)   #Extracting numbers (note error is expected, as r(cde) and r(se_cde) are not being used in my project)
DT_log$var3 <- unlist(str_extract_all(DT_log$txt,  "Imp number|(?<=\\().+?(?=\\))")) #to extract all things within brackets ().Unlist: so output is not a list.
#https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#https://stackoverflow.com/questions/16502219/how-to-use-cast-in-reshape-without-aggregation

DT_log$seq <- with(DT_log, ave(grp, var3, FUN = seq_along))
(gformula_dataset <- dcast(grp + seq ~ var3, data = DT_log, value.var = "var2"))
gformula_dataset$grp <- NULL

#reordering columns
names(gformula_dataset)

col_order <- c("seq", "Imp number", "tce", "se_tce", "nde", "se_nde", "nie", "se_nie", "pm", "se_pm", "cde", "se_cde", "MC_sims", "N")
(void_vol_data <- gformula_dataset[, col_order])


#write.csv(gformula_dataset_1, file = "C:\\Users\\epwkb\\OneDrive - University of Bristol\\MHI\\mediation\\inflammation\\02-analysis\\data\\HPC\\mediation_logs\\UI_ace_score_med_imp.csv", row.names = FALSE) #replace here name as needed.
#===============================================================================
# Histograms of estimates
library(ggplot2)

# TCE
bw <- (2 * IQR(ui_data$tce)) / length(ui_data$tce)^(1/3)

tce <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_histogram( aes(x = tce, y = ..density..), fill="#69b3a2", binwidth = bw ) +
  xlab("TCE") +
  theme_classic()
tce

tce <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = tce, y = ..density..), fill="#69b3a2") +
  xlab("TCE") +
  theme_classic()
tce

se_tce <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = se_tce, y = ..density..), fill="#69b3a2" ) +
  xlab("TCE SE") +
  theme_classic()
se_tce

# NIE
bw <- (2 * IQR(ui_data$nie)) / length(ui_data$nie)^(1/3)

nie <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_histogram( aes(x = nie, y = ..density..), fill="#69b3a2", binwidth = bw ) +
  xlab("NIE") +
  theme_classic()
nie

nie <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = nie, y = ..density..), fill="#69b3a2") +
  xlab("NIE") +
  theme_classic()
nie

se_nie <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = se_nie, y = ..density..), fill="#69b3a2" ) +
  xlab("NIE SE") +
  theme_classic()
se_nie

# NDE
bw <- (2 * IQR(ui_data$nde)) / length(ui_data$nde)^(1/3)

nde <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_histogram( aes(x = nde, y = ..density..), fill="#69b3a2", binwidth = bw ) +
  xlab("NDE") +
  theme_classic()
nde

nde <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = nde, y = ..density..), fill="#69b3a2") +
  xlab("NDE") +
  theme_classic()
nde

se_nde <- ggplot(ui_data, aes(x=x) ) +
  # Top
  geom_density( aes(x = se_nde, y = ..density..), fill="#69b3a2" ) +
  xlab("NDE SE") +
  theme_classic()
se_nde

#===============================================================================
# Set up tibble
# Tibble for estimates
med_mod <- tibble(
  Outcome = character(),
  TCE_coef = numeric(),
  TCE_CI_l = numeric(),
  TCE_CI_u = numeric(),
  NIE_coef = numeric(),
  NIE_CI_l = numeric(),
  NIE_CI_u = numeric(),
  NDE_coef = numeric(),
  NDE_CI_l = numeric(),
  NDE_CI_u = numeric(),
  PM = numeric()
)

#===============================================================================
# Rubins Rules

# UI
med_mod[1,1] <- "Urinary Incontinence"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(ui_data$tce)))
(med_mod[1,2] <- exp(mean(as.numeric(ui_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(ui_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(ui_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(ui_data)+1)/(nrow(ui_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[1,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[1,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(ui_data$nie)))
(med_mod[1,5] <- exp(mean(as.numeric(ui_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(ui_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(ui_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(ui_data)+1)/(nrow(ui_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[1,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[1,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(ui_data$nde)))
(med_mod[1,8] <- exp(mean(as.numeric(ui_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(ui_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(ui_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(ui_data)+1)/(nrow(ui_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[1,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[1,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[1,11] <- mean(as.numeric(ui_data$pm))*100) # mean PM

# Daytime wetting
med_mod[2,1] <- "Daytime wetting"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(dw_data$tce)))
(med_mod[2,2] <- exp(mean(as.numeric(dw_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(dw_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(dw_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(dw_data)+1)/(nrow(dw_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[2,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[2,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(dw_data$nie)))
(med_mod[2,5] <- exp(mean(as.numeric(dw_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(dw_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(dw_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(dw_data)+1)/(nrow(dw_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[2,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[2,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(dw_data$nde)))
(med_mod[2,8] <- exp(mean(as.numeric(dw_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(dw_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(dw_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(dw_data)+1)/(nrow(dw_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[2,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[2,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[2,11] <- mean(as.numeric(dw_data$pm))*100) # mean PM

# Bedtime wetting
med_mod[3,1] <- "Bedtime wetting"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(bw_data$tce)))
(med_mod[3,2] <- exp(mean(as.numeric(bw_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(bw_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(bw_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(bw_data)+1)/(nrow(bw_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[3,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[3,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(bw_data$nie)))
(med_mod[3,5] <- exp(mean(as.numeric(bw_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(bw_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(bw_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(bw_data)+1)/(nrow(bw_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[3,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[3,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(bw_data$nde)))
(med_mod[3,8] <- exp(mean(as.numeric(bw_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(bw_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(bw_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(bw_data)+1)/(nrow(bw_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[3,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[3,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[3,11] <- mean(as.numeric(bw_data$pm))*100) # mean PM

# urgency wetting
med_mod[4,1] <- "Urgency"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(urgency_data$tce)))
(med_mod[4,2] <- exp(mean(as.numeric(urgency_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(urgency_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(urgency_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(urgency_data)+1)/(nrow(urgency_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[4,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[4,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(urgency_data$nie)))
(med_mod[4,5] <- exp(mean(as.numeric(urgency_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(urgency_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(urgency_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(urgency_data)+1)/(nrow(urgency_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[4,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[4,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(urgency_data$nde)))
(med_mod[4,8] <- exp(mean(as.numeric(urgency_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(urgency_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(urgency_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(urgency_data)+1)/(nrow(urgency_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[4,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[4,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[4,11] <- mean(as.numeric(urgency_data$pm))*100) # mean PM

# Nocturia wetting
med_mod[5,1] <- "Nocturia"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(noct_data$tce)))
(med_mod[5,2] <- exp(mean(as.numeric(noct_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(noct_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(noct_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(noct_data)+1)/(nrow(noct_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[5,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[5,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(noct_data$nie)))
(med_mod[5,5] <- exp(mean(as.numeric(noct_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(noct_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(noct_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(noct_data)+1)/(nrow(noct_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[5,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[5,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(noct_data$nde)))
(med_mod[5,8] <- exp(mean(as.numeric(noct_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(noct_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(noct_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(noct_data)+1)/(nrow(noct_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[5,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[5,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[5,11] <- mean(as.numeric(noct_data$pm))*100) # mean PM

# Frequency wetting
med_mod[6,1] <- "Frequency"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(freq_data$tce)))
(med_mod[6,2] <- exp(mean(as.numeric(freq_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(freq_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(freq_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(freq_data)+1)/(nrow(freq_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[6,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[6,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(freq_data$nie)))
(med_mod[6,5] <- exp(mean(as.numeric(freq_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(freq_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(freq_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(freq_data)+1)/(nrow(freq_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[6,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[6,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(freq_data$nde)))
(med_mod[6,8] <- exp(mean(as.numeric(freq_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(freq_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(freq_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(freq_data)+1)/(nrow(freq_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[6,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[6,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[6,11] <- mean(as.numeric(freq_data$pm))*100) # mean PM

# void_post wetting
med_mod[7,1] <- "Voiding postponement"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(void_post_data$tce)))
(med_mod[7,2] <- exp(mean(as.numeric(void_post_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(void_post_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(void_post_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(void_post_data)+1)/(nrow(void_post_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[7,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[7,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(void_post_data$nie)))
(med_mod[7,5] <- exp(mean(as.numeric(void_post_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(void_post_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(void_post_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(void_post_data)+1)/(nrow(void_post_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[7,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[7,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(void_post_data$nde)))
(med_mod[7,8] <- exp(mean(as.numeric(void_post_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(void_post_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(void_post_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(void_post_data)+1)/(nrow(void_post_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[7,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[7,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[7,11] <- mean(as.numeric(void_post_data$pm))*100) # mean PM

# void_vol wetting
med_mod[8,1] <- "Voiding volume"

# TCE [cols 2-4]
(pooledMeanTCE <- mean(as.numeric(void_vol_data$tce)))
(med_mod[8,2] <- exp(mean(as.numeric(void_vol_data$tce)))) # Pooled TCE
(withinVarTCE <- mean(as.numeric(void_vol_data$se_tce^2))) # mean of variances
(betweenVarTCE <- var(as.numeric(void_vol_data$tce))) # variance of coefficients
(dfCorrectionTCE <- (nrow(void_vol_data)+1)/(nrow(void_vol_data))) # dfCorrection
(totVarTCE <- withinVarTCE + betweenVarTCE*dfCorrectionTCE)
#(totVarTCE <- withinVarTCE + betweenVarTCE + betweenVarTCE/80) # same thing
(pooledSETCE <- sqrt(totVarTCE)) # standard error
(med_mod[8,3] <- exp(pooledMeanTCE-(1.96*pooledSETCE))) # lower CI
(med_mod[8,4] <- exp(pooledMeanTCE+(1.96*pooledSETCE))) # upper CI

# NIE [cols 5-7]
(pooledMeanNIE <- mean(as.numeric(void_vol_data$nie)))
(med_mod[8,5] <- exp(mean(as.numeric(void_vol_data$nie)))) # Pooled NIE
(withinVarNIE <- mean(as.numeric(void_vol_data$se_nie^2))) # mean of variances
(betweenVarNIE <- var(as.numeric(void_vol_data$nie))) # variance of coefficients
(dfCorrectionNIE <- (nrow(void_vol_data)+1)/(nrow(void_vol_data))) # dfCorrection
(totVarNIE <- withinVarNIE + betweenVarNIE*dfCorrectionNIE)
#(totVarNIE <- withinVarNIE + betweenVarNIE + betweenVarNIE/80) # same thing
(pooledSENIE <- sqrt(totVarNIE)) # standard error
(med_mod[8,6] <- (exp(pooledMeanNIE)-(1.96*pooledSENIE))) # lower CI
(med_mod[8,7] <- (exp(pooledMeanNIE)+(1.96*pooledSENIE))) # upper CI

# NDE [cols 8-10]
(pooledMeanNDE <- mean(as.numeric(void_vol_data$nde)))
(med_mod[8,8] <- exp(mean(as.numeric(void_vol_data$nde)))) # Pooled NDE
(withinVarNDE <- mean(as.numeric(void_vol_data$se_nde^2))) # mean of variances
(betweenVarNDE <- var(as.numeric(void_vol_data$nde))) # variance of coefficients
(dfCorrectionNDE <- (nrow(void_vol_data)+1)/(nrow(void_vol_data))) # dfCorrection
(totVarNDE <- withinVarNDE + betweenVarNDE*dfCorrectionNDE)
#(totVarNIE <- withinVarNDE + betweenVarNDE + betweenVarNDE/80) # same thing
(pooledSENDE <- sqrt(totVarNDE)) # standard error
(med_mod[8,9] <- exp(pooledMeanNDE-(1.96*pooledSENDE))) # lower CI
(med_mod[8,10] <- exp(pooledMeanNDE+(1.96*pooledSENDE))) # upper CI

# PM [cols 1, 11]
(med_mod[8,11] <- mean(as.numeric(void_vol_data$pm))*100) # mean PM

#===============================================================================
# Make into a pretty table
med_ac_luts <-gt::gt(med_mod, rownames_to_stub = FALSE) %>%
  cols_merge(columns = c("NDE_CI_l", "NDE_CI_u"), pattern = "{1}, {2}") %>%
  cols_merge(columns = c("NIE_CI_l", "NIE_CI_u"), pattern = "{1}, {2}") %>%
  cols_merge(columns = c("TCE_CI_l", "TCE_CI_u"), pattern = "{1}, {2}") %>%
  fmt_number(decimals = 2) %>%
  cols_align(align = "center") %>%
  cols_align(align = "left", columns = "Outcome") %>%
  cols_label(
    NDE_coef = "OR",
    NDE_CI_l = "95% CI",
    NIE_coef = "OR",
    NIE_CI_l = "95% CI",
    TCE_coef = "OR",
    TCE_CI_l = "95% CI"
  ) %>%
  tab_spanner(
    label = "Natural direct effect",
    columns = c("NDE_coef", "NDE_CI_l")
  ) %>%
  tab_spanner(
    label = "Natural indirect effect",
    columns = c("NIE_coef", "NIE_CI_l")
  ) %>%
  tab_spanner(
    label = "Total causal effect",
    columns = c("TCE_coef", "TCE_CI_l")
  )


med_ac_luts

gtsave(med_ac_luts, file = paste0(Tables,"table_3.docx"))

# END
