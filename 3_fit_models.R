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

#dataset without high PM outliers
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95)
cco_low_PM <- cco[pm25_3day < pm25_3day_q95]

# path for model results
model_path <- paste0("results/models/")


#---------- set up for model fitting ----------#

# # run this for parallel computing (all chains can run at the same time)
# options(mc.cores = parallel::detectCores())

# run this to run one chain at a time
options(mc.cores = 1)

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

