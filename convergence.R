###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Check model convergence
# Inputs: 
#   Model output
# Output:
#   Convergence plots (r-hat values)
###############################################################

library(bayesplot)
library(rstan)
library(rstanarm)
library(brms)
library(ggplot2)


# path to store convergence results
converge_path <- "data/intermediate/convergence/"

# model names
mod.names <- c("spline_trim95",
               "t2_trim95",
               "spline_trim99",
               "spline_trim95_temp3day")


# loop through models
for (i in 1:length(mod.names)) {
  
  cat("Running: ", mod.names[i], "\n")
  
  #----- load model ------#
  
  cat("Load model\n")
  load(paste0("results/models/", mod.names[i], ".RData"))
  gc()
  
  #----- get rhats -----#
  
  # (takes a long time to run, so save results)
  cat("Get rhats\n")
  rhats <- brms::rhat(stan.mod)
  save(rhats, file = paste0(converge_path, "rhats_", mod.names[i], ".RData"))
  
  
  #----- plot -----#
  
  cat("Plot rhats\n")
  bayesplot::mcmc_rhat(rhats)
  ggsave(paste0(converge_path, "rhat_plot_", mod.names[i], ".pdf"), 
         width = 8, height = 6)
  
  gc()
}
