###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Extract results from a single STAN model
# Input: 
#   Output from one model
# Outputs: 
#   Model results (figures and tables)
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
library(brms)
library(cowplot)
library(patchwork)


#---------- user input ----------#

### choose model name (run this script once for each model to get results)

mod.name <- "spline_trim95"
#mod.name <- "t2_trim95"
#mod.name <- "spline_trim99"
#mod.name <- "spline_trim95_temp3day"


#---------- path for figures ----------#

# path for saving figures
fig_path <- paste0("results/figures/tmax_temp_", mod.name, "/")

# create directory if it doesn't exist yet
if(mod.name == "spline_trim95" & !dir.exists(fig_path)){
  dir.create(fig_path)
}

#---------- load case-crossover data ----------#

# Load case-crossover data (main analysis)

load(paste0("data/intermediate/cco_trim95_tmax_temp.RData"))
cco <- cco_trim95
rm(cco_trim95)


#---------- load model results ----------#

# load model
load(paste0("results/models/", mod.name, ".RData"))
gc()


#---------- define counterfactuals of interest for RERI ----------#

# these contrasts are the same for all models for direct comparison

# get median & 95th percentile of exposures for main analysis
med_temp <- median(cco[cases == 1, temp_lag0])
q95_temp <- quantile(cco[cases == 1, temp_lag0], 0.95)
med_pm25 <- median(cco[cases == 1, pm25_3day])
q95_pm25 <- quantile(cco[cases == 1, pm25_3day], 0.95)

# set contrasts (unexposed and exposed) for temperature and PM2.5
temp_exposed <- q95_temp
temp_unexposed <- med_temp
pm25_exposed <- q95_pm25
pm25_unexposed <- med_pm25


#---------- if trimming at 99th percentile, load that dataset instead ----------#

if (mod.name == "spline_trim99") {
  load(paste0("data/intermediate/cco_trim99_tmax_temp.RData"))
  cco <- cco_trim99
  rm(cco_trim99)
}


#---------- get basis/X values to use for counterfactual exposures ----------#

# write a function to get cf11, cf10, cf01, cf00 matrices (spline basis)
get_cf_basis <- function(temp, pm25){
  
  # use model fit to get design matrix for new data
  
  # use same-day or three-day temperature (depends on model)
  if (mod.name == "spline_trim95_temp3day") {
    new_expo <- data.frame(temp_3day = temp, 
                           pm25_3day = pm25,
                           cases = 1, 
                           keyid = first(cco$keyid)) # can be any ID
  } else {
    new_expo <- data.frame(temp_lag0 = temp, 
                           pm25_3day = pm25,
                           cases = 1, 
                           keyid = first(cco$keyid)) # can be any ID
  }
  
  newX <- brms:::standata.brmsfit(stan.mod, 
                                  newdata = new_expo,
                                  allow_new_levels = TRUE)
  design.mat <- cbind(newX$X, newX$Xs)
  
  return(design.mat)
  
}

# now apply function to get each basis
cf11 <- get_cf_basis(temp = temp_exposed, pm25 = pm25_exposed)
cf10 <- get_cf_basis(temp = temp_exposed, pm25 = pm25_unexposed)
cf01 <- get_cf_basis(temp = temp_unexposed, pm25 = pm25_exposed)
cf00 <- get_cf_basis(temp = temp_unexposed, pm25 = pm25_unexposed)


#---------- get posterior betas ----------#

# get model betas (posterior distribution)
beta.post <- as_draws_df(stan.mod)

# select columns that start with "b" (ignore sd and random effects)
beta.post <- beta.post |>
  select(starts_with("b")) |>
  as.matrix() |>
  suppressWarnings()
  

#---------- estimate ORs & RERI for each counterfactual ----------#

# initialize data frame
mod.res <- as.data.frame(matrix(nrow = nrow(beta.post), ncol = 4)) |>
  setNames(c("OR11", "OR10", "OR01", "RERI"))

# loop through rows (draws from posterior distribution)
for(i in 1:nrow(beta.post)){
  
  # get odds ratio comparing 11 to reference
  mod.res$OR11[i] <- exp(cf11 %*% beta.post[i,] - cf00 %*% beta.post[i,])
  
  # get odds ratio comparing 10 to reference
  mod.res$OR10[i] <- exp(cf10 %*% beta.post[i,] - cf00 %*% beta.post[i,])
  
  # get odds ratio comparing 01 to reference
  mod.res$OR01[i] <- exp(cf01 %*% beta.post[i,] - cf00 %*% beta.post[i,])
  
  # calculate RERI
  mod.res$RERI[i] <- mod.res$OR11[i] - mod.res$OR10[i] - mod.res$OR01[i] + 1
  
}

# now just take column means/quantiles (across draws)
point.est <- apply(mod.res, 2, function(x) mean(x))
ci95 <- t(apply(mod.res, 2, function(x) quantile(x, c(0.025, 0.975))))

# get posterior probabilities
prob <- c(apply(mod.res[,1:3], 2, function(x) mean(x <= 1)), # null is 1 for ORs
           mean(mod.res[,4] < 0)) # null is 0 for RERI

# set up data for ggplot
res_df <- cbind.data.frame(point.est, ci95, prob)

# column with more names for plotting (in an order that makes sense)
res_df$estimand <- rownames(res_df)
res_df$estimand <- factor(res_df$estimand, 
                                levels = rev(c("OR10", 
                                               "OR01", 
                                               "OR11", 
                                               "RERI")))

res_df <- res_df %>%
  arrange(desc(estimand))

# create a column for the model name
res_df$model <- mod.name

# save this table (results saved here will be merged across models in script 5)
saveRDS(res_df, file = paste0("data/intermediate/OR_tables/OR_table_tmax_temp_", 
                              mod.name, ".rds"))


#---------- get interaction on multiplicative scale ----------#

# get the interaction on the multiplicative scale (exponentiate the interaction beta)
if(str_detect(mod.name, "spline")){
  
  # name of interaction coefficient depends on model
  if (mod.name == "spline_trim95_temp3day"){
    mult <- beta.post[,"b_temp_3day:pm25_3day"]
  } else {
    mult <- beta.post[,"b_temp_lag0:pm25_3day"]
  }
  
  mult_mean <- exp(mean(mult))
  mult_ci <- exp(quantile(mult, c(0.025, 0.975)))
  paste0(mult_mean, " (95% CI: ", mult_ci[1], ", ", mult_ci[2], ")")
}

#---------- visualize effect of just temp ----------#

# how does OR10 change for different temps when we hold PM at the median?

# reference temp
ref_temp <- med_temp

# vector of comparison temps
comp_temp_vec <- seq(60, 115, by = 1)

# constant PM2.5
constant_pm25 <- med_pm25

# initialize data frame for results
temp_res <- as.data.frame(matrix(nrow = length(comp_temp_vec), ncol = 4)) |>
  setNames(c("comp_expo", "OR", "lower", "upper"))
temp_res$comp_expo <- comp_temp_vec

# loop through comparison temperatures
for(j in 1:length(comp_temp_vec)){
  
  comp_expo <- comp_temp_vec[j]
  
  # now apply function to get each counterfactual
  cf10_temp <- get_cf_basis(temp = comp_expo, pm25 = constant_pm25)
  cf00_temp <- get_cf_basis(temp = ref_temp, pm25 = constant_pm25)
  
  # initialize data frame
  OR <- vector()
  
  # loop through rows and do the calc above
  for(i in 1:nrow(beta.post)){
    
    # get odds ratio comparing 10 to reference
    OR[i] <- exp(cf10_temp %*% beta.post[i,] - cf00_temp %*% beta.post[i,])
    
  }
  
  # now just take column means/quantiles
  point.est <- mean(OR)
  ci95 <- quantile(OR, c(0.025, 0.975))
  
  temp_res[j,2:4] <- c(point.est, ci95)
  
}

# get temperature in celsius
temp_res <- temp_res %>%
  mutate(comp_expo = (comp_expo - 32) * (5/9))


#---------- visualize effect of just PM2.5 ----------#

# reference PM2.5
ref_pm25 <- med_pm25

# vector of comparison PM2.5 values
comp_pm25_vec <- seq(1, 18, by = 1)

# constant temp
constant_temp <- med_temp

# initialize data frame for results
pm25_res <- as.data.frame(matrix(nrow = length(comp_pm25_vec), ncol = 4)) |>
  setNames(c("comp_expo", "OR", "lower", "upper"))
pm25_res$comp_expo <- comp_pm25_vec

# loop through comparison PM2.5 values
for(j in 1:length(comp_pm25_vec)){
  
  comp_expo <- comp_pm25_vec[j]
  
  # now apply function to get each counterfactual
  cf01 <- get_cf_basis(temp = constant_temp, pm25 = comp_expo)
  cf00 <- get_cf_basis(temp = constant_temp, pm25 = ref_pm25)
  
  # initialize data frame
  OR <- vector()
  
  # loop through rows and do the calc above
  for(i in 1:nrow(beta.post)){
    
    # get odds ratio comparing 10 to reference
    OR[i] <- exp(cf01 %*% beta.post[i,] - cf00 %*% beta.post[i,])
    
  }
  
  # now just take column means/quantiles
  point.est <- mean(OR)
  ci95 <- quantile(OR, c(0.025, 0.975))
  
  pm25_res[j,2:4] <- c(point.est, ci95)
  
}


#---------- plot temp and PM together ----------#

# if main results, save plot
if (mod.name == "spline_trim95") {
  
  # temp plot
  temp_effect_plot <- temp_res |>
    ggplot(aes(x = comp_expo, y = OR)) +
    geom_hline(yintercept = 1, col = "black", lty = 1, linewidth = 0.3) +
    geom_line(col = "#ffad4d") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#ffad4d") +  
    geom_vline(xintercept = (med_temp - 32) * (5/9), lty = 2, col = "gray40") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.82, 1.25)) +
    labs(x = "Daily maximum temperature (\u00B0C)",
         y = "Odds ratio (relative to median)") +
    theme_minimal() +
    theme(panel.grid.major = element_line(linewidth = 0.2))
  
  # pm25 plot
  pm25_effect_plot <- pm25_res |>
    ggplot(aes(x = comp_expo, y = OR)) +
    geom_hline(yintercept = 1, col = "black", lty = 1, linewidth = 0.3) +
    geom_line(col = "#631ea2") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#631ea2") +
    geom_vline(xintercept = med_pm25, lty = 2, col = "gray40") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(#breaks = c(0.95, 1, 1.05),
      expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.84, 1.25)) +
    labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
         y = "Odds ratio (relative to median)") +
    theme_minimal() +
    theme(panel.grid.major = element_line(linewidth = 0.2))
  
  # align plots (so the plotting space is the same)
  temp_pm25_plots <- cowplot::align_plots(temp_effect_plot, pm25_effect_plot, align = "hv")
  
  # plot
  plot_grid(temp_pm25_plots[[1]], temp_pm25_plots[[2]])
  ggsave(paste0(fig_path, "independent_effects.pdf"), height = 3, width = 6.5)
  
# otherwise save plotting data. all sensitivity analysis plotted together in script 5
} else {
  
  # save plotting data
  temp_res$model <- mod.name
  pm25_res$model <- mod.name
  
  saveRDS(temp_res, 
          file = paste0("data/intermediate/independent_effects/independent_effects_temp_tmax_temp_", 
                        mod.name, ".rds"))
  
  saveRDS(pm25_res, 
          file = paste0("data/intermediate/independent_effects/independent_effects_pm25_tmax_temp_", 
                        mod.name, ".rds"))
  
}


#---------- risk surface plot ----------#

# get pairs of PM and temperature values 
temp_vals <- seq(60, 115, length = 10)
pm25_vals <- seq(1, 18, length = 10)
temp_pm25_grid <- expand.grid(temp = temp_vals, pm25 = pm25_vals)
temp_pm25_grid$OR11 <- NA

# loop through temperature and PM2.5 pairs
for(i in 1:nrow(temp_pm25_grid)){
  
  # get current temp and PM2.5
  temp_i <- temp_pm25_grid$temp[i]
  pm25_i <- temp_pm25_grid$pm25[i]
  
  # get basis for current values
  cf11 <- get_cf_basis(temp = temp_i, pm25 = pm25_i)
  cf00 <- get_cf_basis(temp = med_temp, pm25 = med_pm25)
  
  # initialize vector for odds ratios
  # comparing both exposures at current values to both exposures at median
  OR11 <- c()
  
  # loop through beta.post to get OR11 estimates
  for(j in 1:nrow(beta.post)){
    # get odds ratio comparing 11 to reference
    OR11[j] <- exp(cf11 %*% beta.post[j,] - cf00 %*% beta.post[j,])
  }
  
  # get posterior mean
  temp_pm25_grid$OR11[i] <- mean(OR11)
    
}

# get temperature in celsius
temp_pm25_grid <- temp_pm25_grid %>%
  mutate(temp = (temp - 32) * (5/9))

# color scale
fill_scale <- scale_fill_gradient2(low = "#056608", mid = "white", high = "#a73672",
                                   midpoint = 1,
                                   limits = c(0.84, 1.331)
                                   )

# for main results, save plot
if (mod.name == "spline_trim95") {
  
  temp_pm25_grid |>
    ggplot() + 
    geom_tile(aes(x = temp, y = pm25, fill = OR11)) +
    labs(x = "Daily maximum temperature (\u00B0C)",
         y = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
         fill = "Odds ratio (relative to \u00d7)"
    ) +
    fill_scale +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_point(x = 29.6, y = 8.9, shape = 4, size = 2, stroke = 1) +
    theme_minimal() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_colorbar(title.position = "top",
                                 barwidth = 15,
                                 barheight = 0.5,
                                 ticks.colour = NA))
  ggsave(paste0(fig_path, "risk_surface.pdf"), width = 4, height = 4.5)
  
  
  
# otherwise save data and plot all sensitivity analysis results together in script 5
} else {
  
  # save plotting data
  temp_pm25_grid$model <- mod.name
  
  saveRDS(temp_pm25_grid, 
          file = paste0("data/intermediate/risk_surface/risk_surface_tmax_temp_", 
                        mod.name, ".rds"))
  
}

