###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Fit models with STAN
# Input: 
#   Case-crossover data
# Outputs: 
#   Models
###############################################################

# Load libraries
library(data.table)
library(survival)
library(rstanarm)
library(shinystan)
library(ggplot2)
library(viridis)
library(ggpubr)
library(tidyverse)
library(splines)
library(mgcv)
library(brms)

#---------- user inputs ----------#

# choose temperature metric
dataset <- "tmax_temp"
#dataset <- "tmax_pcnt"
#dataset <- "hi_temp"
#dataset <- "hi_pcnt"


#---------- load case-crossover data ----------#

# Load trimmed case-crossover data (not trimmed)
# (this way I can always change how we trim for sensitivity analyses later)
load(paste0("data/intermediate/cco_", dataset, ".RData"))

# Sort data by strata (required for STAN)
cco <- cco[order(as.factor(cco$bene_id)),]

# dataset without high PM outliers
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95)
cco_low_PM <- cco[pm25_3day < pm25_3day_q95]

# path for model results
model_path <- paste0("results/models/")


#---------- set up for model fitting ----------#

# run this for parallel computing (all chains can run at the same time)
options(mc.cores = parallel::detectCores())

# # run this to run one chain at a time (uses less memory but takes longer)
# options(mc.cores = 1)

# function to fit & save models with brms
fit_brms_model <- function(mod.data,
                           mod.form,
                           mod.name){
  
  # get file name to save model results
  mod.file <- paste0(model_path, mod.name, ".RData")
  
  # fit model
  stan.mod <- brm(mod.form, 
                  data = mod.data,
                  family = bernoulli(link = "logit"))
  save(stan.mod, file = mod.file)
}

# function to fit & save models with rstanarm
fit_rstanarm_model <- function(mod.data,
                               mod.form,
                               mod.name){
  
  # get file name to save model results
  mod.file <- paste0(model_path, mod.name, ".RData")
  
  # fit model
  stan.mod <- stan_clogit(mod.form, 
                          strata = bene_id, 
                          data = mod.data)
  save(stan.mod, file = mod.file)
}

#---------- model fitting ----------#

# main analysis
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) +
  temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "spline_lowPM")
gc()

# sensitivity analysis
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ t2(temp_lag0, pm25_3day) + (1|bene_id)"),
  mod.name = "t2_lowPM")
gc()







#################################################################
###################### OLD MODELS ###############################
#################################################################


#------------------------- newest models with brms

# s (instead of ns) for main effects, linear interaction
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ s(temp_lag0) + s(pm25_3day) + temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "NEW_brms_s_linear")
gc()

# use s AND t2
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ s(temp_lag0) + s(pm25_3day) + t2(temp_lag0, pm25_3day) + (1|bene_id)",
  mod.name = "NEW_brms_s_t2")
gc()

# use ns AND t2
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) + t2(temp_lag0, pm25_3day) + (1|bene_id)",
  mod.name = "NEW_brms_ns_t2")
gc()




#------------------------- models with rstanarm

##### MAIN ANALYSIS #####

# rstanarm, full dataset, linear interaction
fit_rstanarm_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) + 
  temp_lag0:pm25_3day"),
  mod.name = "rstanarm_spline_lowPM_all")

# this model finished running in nine days


##### SENSITIVITY ANALYSIS #####

# t2 tensor product isn't compatible with rstanarm
# using code from Federica to create basis manually

# create basis manually
joint_basis <- function(x1, x2, df = 3){
  b1 <- bs(x1, df=df)
  b2 <- bs(x2, df=df)
  basis <- cbind(b1, b2,
                 b1[,1]*b2[,1], b1[,1]*b2[,2], b1[,1]*b2[,3],
                 b1[,2]*b2[,1], b1[,2]*b2[,2], b1[,2]*b2[,3],
                 b1[,3]*b2[,1], b1[,3]*b2[,2], b1[,3]*b2[,3])
  return(list(b1 = b1, b2 = b2, basis = basis))
}

# use the functions above to get the basis as columns in the data
basis_el <- joint_basis(cco_low_PM$temp_lag0, cco_low_PM$pm25_3day)
b1 <- basis_el$b1
b2 <- basis_el$b2
basis_matrix <- basis_el$basis
colnames(basis_matrix) <- c(paste0("b1_", 1:3), 
                            paste0("b2_", 1:3), 
                            paste0("int_", 1:9))
cco_low_PM <- cbind(cco_low_PM, basis_matrix)

# rstanarm, full dataset, non-linear interaction
fit_rstanarm_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ b1_1 + b1_2 + b1_3 + b2_1 + b2_2 + b2_3 +
                        int_1 + int_2 + int_3 + int_4 + int_5 + 
                        int_6 + int_7 + int_8 + int_9"), 
  mod.name = "rstanarm_manual_basis_lowPM_all")

#### ^^ this model will take months to finish running (more covariates than main model)



# try brms with te and t2 as Kevin suggested
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ s(temp_lag0) + s(pm25_3day) + 
                        t2(temp_lag0, pm25_3day) + (1|bene_id)"),
  mod.name = "brms_te_t2")





#################################################################
###################### OLD(er) MODELS ###########################
#################################################################



#---------- brms ----------#

# # full dataset, exclude high PM, linear interaction
# fit_brms_model(
#   mod.data = cco_low_PM,
#   mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) + 
#   temp_lag0:pm25_3day + (1|bene_id)",
#   mod.name = "spline_lowPM")
# gc()
# 
# # # full dataset, exclude high PM, tensor with s
# # fit_brms_model(
# #   mod.data = cco_low_PM,
# #   mod.form = as.formula("cases ~ s(temp_lag0, pm25_3day) + (1|bene_id)"),
# #   mod.name = "s_lowPM")
# 
# # try with t2 instead of s!
# # full dataset, exclude high PM, tensor with t2
# fit_brms_model(
#   mod.data = cco_low_PM,
#   mod.form = as.formula("cases ~ t2(temp_lag0, pm25_3day) + (1|bene_id)"),
#   mod.name = "t2_lowPM")


#---------- compare models with rstanarm and brms (data subset) ----------#

# use random sample of data for noW!!
set.seed(17)
id_keep <- sample(unique(cco_low_PM$bene_id), 30000)
cco_low_PM <- cco_low_PM %>%
  filter(bene_id %in% id_keep)

# note 5 df
fit_rstanarm_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ ns(temp_lag0, df = 5) + ns(pm25_3day, df = 5) + 
  temp_lag0:pm25_3day"),
  mod.name = "rstanarm_spline_lowPM_seed17")
gc()

# and brms for comparison
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 5) + ns(pm25_3day, df = 5) +
  temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "brms_spline_lowPM_seed17")
gc()

# # I can't use rstanarm and t2, but can I use s?
# fit_rstanarm_model(
#   mod.data = cco_low_PM,
#   mod.form = as.formula("cases ~ s(temp_lag0, pm25_3day)"),
#   mod.name = "rstanarm_s_lowPM")
# # no.


# what if we only remove the top 1% of PM2.5?

pm25_3day_q99 <- quantile(cco$pm25_3day, 0.99)
cco_low_PM99 <- cco[pm25_3day < pm25_3day_q99]

set.seed(17)
id_keep99 <- sample(unique(cco_low_PM99$bene_id), 30000)
cco_low_PM99 <- cco_low_PM99 %>%
  filter(bene_id %in% id_keep99)

# bottom 99% of PM
fit_rstanarm_model(
  mod.data = cco_low_PM99,
  mod.form = as.formula("cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) + 
  temp_lag0:pm25_3day"),
  mod.name = "rstanarm_spline_lowPM99")
gc()

# more df
fit_rstanarm_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ ns(temp_lag0, df = 5) + ns(pm25_3day, df = 5) + 
  temp_lag0:pm25_3day"),
  mod.name = "rstanarm_spline_lowPM_5df")


################################################################

# other experimenting


# spline for dataset with raw temp or HI
model_form_spline_temp <- "cases ~ ns(temp_lag0, knots = c(80, 90)) +
ns(pm25_3day, knots = c(7, 12)) + temp_lag0:pm25_3day"

# full dataset, exclude high PM, linear interaction
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, knots = c(80, 90)) +
ns(pm25_3day, knots = c(7, 12)) + temp_lag0:pm25_3day",
  mod.name = "spline_lowPM")
gc()

# full dataset, exclude high PM, tensor with t2
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ t2(temp_lag0, pm25_3day) + (1|bene_id)"),
  mod.name = "t2_lowPM")

# full dataset, exclude high PM, linear interaction
fit_brms_model(
  mod.data = cco,
  mod.form = "cases ~ ns(temp_lag0, knots = c(80, 90)) +
ns(pm25_3day, knots = c(7, 12)) + temp_lag0:pm25_3day",
  mod.name = "spline_allPM")
gc()

# full dataset, all PM, linear interaction with splines but more knots
fit_brms_model(
  mod.data = cco,
  mod.form = as.formula("cases ~ ns(temp_lag0, knots = c(70, 80, 90, 95, 100)) +
ns(pm25_3day, knots = c(5, 8, 15, 25, 75)) + temp_lag0:pm25_3day"),
  mod.name = "spline_allPM_5knots")
gc()

# full dataset, all PM, linear interaction with splines but more knots
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ ns(temp_lag0, knots = c(70, 80, 90, 95, 105)) +
ns(pm25_3day, knots = c(5, 8, 11, 13, 16)) + temp_lag0:pm25_3day"),
  mod.name = "spline_lowPM_5knots")
gc()


#########################################################################################

# some experimenting

# full dataset, exclude high PM, linear interaction, NEW KNOTS, same-day PM
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_lag0, df = 3) + temp_lag0:pm25_lag0",
  mod.name = "spline_lowPM_newKnots_lag0")
gc()

# full dataset, exclude high PM, linear interaction, more df
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 4) + ns(pm25_3day, df = 4) + temp_lag0:pm25_3day",
  mod.name = "spline_lowPM_4df")
gc()

#########################################################################################

#----- fit final model on subgroups

# medication use 
antichol <- cco_low_PM[user_Antichol == 1]
anticholX <- cco_low_PM[is.na(user_Antichol) | user_Antichol == 0]
stim <- cco_low_PM[user_Stim == 1]
stimX <- cco_low_PM[is.na(user_Stim) | user_Stim == 0]
loopd <- cco_low_PM[user_LoopD == 1]
loopdX <- cco_low_PM[is.na(user_LoopD) | user_LoopD == 0]

# Medicaid eligibility
dual <- cco_low_PM[dualeligible == 1]
dualX <- cco_low_PM[dualeligible == 0]

# age
old <- cco_low_PM[age > median(cco_low_PM$age)]
oldX <- cco_low_PM[age <= median(cco_low_PM$age)]

# sex
male <- cco_low_PM[sex == 1]
maleX <- cco_low_PM[sex == 2]

# race
rWhite <- cco_low_PM[race == 1]
rBlack <- cco_low_PM[race == 2]
#rOther <- cco[race == 3]
rAsian <- cco_low_PM[race == 4]
rHispanic <- cco_low_PM[race == 5]
#rNative <- cco[race == 6]



#----- drugs

fit_brms_model(
  mod.data = antichol,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_antichol")
gc()

fit_brms_model(
  mod.data = anticholX,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_anticholX")
gc()


fit_brms_model(
  mod.data = stim,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_stim")
gc()

fit_brms_model(
  mod.data = stimX,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_stimX")
gc()

fit_brms_model(
  mod.data = loopd,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_loopd")
gc()

fit_brms_model(
  mod.data = loopdX,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_loopdX")
gc()


#----- race

fit_brms_model(
  mod.data = rWhite,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_rWhite")
gc()

fit_brms_model(
  mod.data = rBlack,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_rBlack")
gc()

fit_brms_model(
  mod.data = rAsian,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_rAsian")
gc()

fit_brms_model(
  mod.data = rHispanic,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_rHispanic")
gc()


#----- age

fit_brms_model(
  mod.data = old,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_old")
gc()

fit_brms_model(
  mod.data = oldX,
  mod.form = model_form_spline_temp,
  mod.name = "spline_lowPM_strat_oldX")
gc()
