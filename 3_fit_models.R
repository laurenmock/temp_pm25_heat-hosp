###############################################################
# Heat-related outcomes
# Author: Lauren Mock
# Aim: Fit models with STAN
# Input: Case-crossover data
# Outputs: Models
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

# Load case-crossover data (not trimmed)
load(paste0("data/intermediate/cco_", dataset, ".RData"))

# Sort data by strata (required for STAN)
cco <- cco[order(as.factor(cco$bene_id)),]

# mini test dataset
cco_mini <- head(cco, 100)

#dataset without high PM outliers
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95)
cco_low_PM <- cco[pm25_3day < pm25_3day_q95]

# path for model results
model_path <- paste0("results/models/")



###############################################################################

#---------- function to fit models with brms & save ----------#

# need to run this for parallel computing (all chains can run at the same time)
options(mc.cores = parallel::detectCores())

# function to fit & save models
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


#---------- fit models with brms ----------#

##### FINAL MODELS #####

# full dataset, exclude high PM, linear interaction
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) + 
  temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "spline_lowPM")
gc()

# full dataset, exclude high PM, tensor with s
fit_brms_model(
  mod.data = cco_low_PM,
  mod.form = as.formula("cases ~ s(temp_lag0, pm25_3day) + (1|bene_id)"),
  mod.name = "s_lowPM")


################################################################

# other experimenting
### FEDERICA, don't worry about anything below this!

################################################################


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
