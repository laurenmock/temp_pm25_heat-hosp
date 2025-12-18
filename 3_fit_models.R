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
library(tidyverse)
library(splines)
library(mgcv)
library(brms)


#---------- load case-crossover data ----------#

# Load case-crossover datasets
load(paste0("data/intermediate/cco_trim95_tmax_temp.RData"))
load(paste0("data/intermediate/cco_trim99_tmax_temp.RData"))

# Sort data by strata (required for STAN)
cco_trim95 <- cco_trim95[order(as.factor(cco_trim95$bene_id)),]
cco_trim99 <- cco_trim99[order(as.factor(cco_trim99$bene_id)),]

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


#---------- MODEL FITTING ----------#

# main analysis
fit_brms_model(
  mod.data = cco_trim95,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) +
  temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "spline_trim95")
gc()

# sensitivity analysis with tensor product
fit_brms_model(
  mod.data = cco_trim95,
  mod.form = as.formula("cases ~ t2(temp_lag0, pm25_3day) + (1|bene_id)"),
  mod.name = "t2_trim95")
gc()

# sensitivity analysis with trimming at 99th percentile
fit_brms_model(
  mod.data = cco_trim99,
  mod.form = "cases ~ ns(temp_lag0, df = 3) + ns(pm25_3day, df = 3) +
  temp_lag0:pm25_3day + (1|bene_id)",
  mod.name = "spline_trim99")
gc()

# sensitivity analysis with three-day temperature exposure
fit_brms_model(
  mod.data = cco_trim95,
  mod.form = "cases ~ ns(temp_3day, df = 3) + ns(pm25_3day, df = 3) +
  temp_3day:pm25_3day + (1|bene_id)",
  mod.name = "spline_trim95_temp3day")
gc()

