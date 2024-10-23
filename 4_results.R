###############################################################
# Heat-related outcomes
# Author: Lauren Mock
# Aim: Extract results from STAN models (spline effects)
# Input: Models
# Outputs: Model results
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


#---------- user inputs ----------#

# choose temperature metric
dataset <- "tmax_temp"
#dataset <- "tmax_pcnt"
#dataset <- "hi_temp"
#dataset <- "hi_pcnt"

### choose model name

#--- splines with linear interaction
mod.name <- "spline_lowPM.RData"
#mod.name <- "spline_lowPM_newKnots.RData"


#--- tensor product (complex interaction)
#mod.name <- "s_lowPM.RData"


#--- drugs
#mod.name <- "spline_lowPM_strat_anticholX.RData"
#mod.name <- "spline_lowPM_strat_antichol.RData"
#mod.name <- "spline_lowPM_strat_loopdX.RData"
#mod.name <- "spline_lowPM_strat_loopd.RData"

#--- race
#mod.name <- "spline_lowPM_strat_rAsian.RData"

#--- age
#mod.name <- "spline_lowPM_strat_oldX.RData"



#--- other things
#mod.name <- "spline_lowPM_newKnots_lag0.RData"
#mod.name <- "spline_lowPM_4df.RData"


#---------- path for figures ----------#

# path for saving figures
fig_path <- paste0("results/figures/", dataset, "_", 
                   str_remove(mod.name, ".RData"), "/")

# create directory if it doesn't exist yet
if(!dir.exists(fig_path)){
  dir.create(fig_path)
}

#---------- load case-crossover data ----------#

###### NOTE: 
# I'm using the dataset before trimming to get median, percentiles, etc.
# for both allPM and low PM models. 

# Load case-crossover data
load(paste0("data/intermediate/cco_trimmed_", dataset, ".RData"))


#---------- load model results ----------#

# path where model results are stored
mod.path <- paste0("results/models/")

# load model
load(paste0(mod.path, mod.name))
gc()

# use something like this if loading multiple files at once 
#spline.files <- list.files(model_path, pattern = "*spline")
# spline.files <- list.files(model_path,
#                            pattern = dataset
# )
# 
# spline.files <- spline.files[str_detect(spline.files, "spline")]

# spline.files <- paste0(model_path, file)
# mod.name <- file

#---------- define counterfactuals of interest for RERI ----------#

# get median & 95th percentile of exposures
med_temp <- median(cco[cases == 1, temp_lag0])
q95_temp <- quantile(cco[cases == 1, temp_lag0], 0.90)
med_pm25 <- median(cco[cases == 1, pm25_lag0])
q95_pm25 <- quantile(cco[cases == 1, pm25_3day], 0.90)

# set "unexposed" and "exposed" levels for each exposure
temp_exposed <- q95_temp
temp_unexposed <- med_temp
pm25_exposed <- q95_pm25
pm25_unexposed <- med_pm25


#---------- get basis/X values to use for counterfactual exposures ----------#

# write a function to get cf11, cf10, cf01, cf00 matrices (spline basis)
# plug in temp and PM2.5 values
get_cf_basis <- function(temp, pm25){
  
  # use model fit to get design matrix for new data
  new_expo <- data.frame(temp_lag0 = temp, 
                         pm25_3day = pm25,
                         cases = 1, 
                         keyid = first(cco$keyid))
  
  # for testing PM2.5 lag0
  if(str_detect(mod.name, "lag0")){
    new_expo <- data.frame(temp_lag0 = temp, pm25_lag0 = pm25,
                           cases = 1, keyid = first(cco$keyid))
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

# loop through rows (draws from posterior distribution) and do the calc above
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

# get p-values
pval <- c(apply(mod.res[,1:3], 2, function(x) mean(x <= 1)),
           mean(mod.res[,4] < 0))

# set up data for ggplot
res_df <- cbind.data.frame(point.est, ci95, pval)

# column with more descriptive names for plotting (in an order that makes sense)
res_df$estimand <- c("temp_pm25", 
                     "temp", 
                     "pm25", 
                     "RERI")
res_df$estimand <- factor(res_df$estimand, 
                                levels = rev(c("temp", 
                                               "pm25", 
                                               "temp_pm25", 
                                               "RERI")))

# column for facet
res_df$groups <- c("Odds ratios", "Odds ratios", 
                         "Odds ratios", "Interaction")

res_df$groups <- factor(res_df$groups, levels = c("Odds ratios", "Interaction"))

# column for vertical line at null (1 for OR, 0 for RERI)
res_df$null.val <- c(1,1,1,0)

# # or plot with both (facet_wrap)
# res_df |>
#   ggplot(aes(x = point.est, y = estimand)) +
#   geom_point() +
#   geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0) +
#   geom_vline(aes(xintercept = null.val), col = "gray60", lty = 2) +
#   labs(title = "", x = "", y = "") +
#   facet_wrap(~groups, nrow = 1, scales = "free") +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text = element_text(size = 11),
#         plot.margin = margin(1, 15, 1, 1))

or_plot <- res_df |>
  filter(groups == "Odds ratios") %>%
  ggplot(aes(x = point.est, y = estimand)) +
  geom_vline(aes(xintercept = null.val), col = "gray60", lty = 2) +
  geom_point() +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  labs(title = "Odds ratios", x = "", y = "") +
  xlim(c(0.99, 1.11)) +
  #scale_y_discrete(labels = function(X) parse(text = x)) +
  scale_y_discrete(labels = c("temp" = expression("Exposure to heat (" * OR[10] * ")"), 
                              "pm25" = expression("Exposure to " * PM[2.5] * " (" * OR["01"] * ")"), 
                              "temp_pm25" = expression("Exposure to heat & " 
                                                       * PM[2.5] * " (" * OR[11] * ")"))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.margin = margin(1, 15, 1, 1))

reri_plot <- res_df |>
  filter(groups == "Interaction") %>%
  ggplot(aes(x = point.est, y = estimand)) +
  geom_vline(aes(xintercept = null.val), col = "gray60", lty = 2) +
  geom_point() +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  labs(title = "Interaction", x = "", y = "") +
  xlim(-0.001, 0.036) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.margin = margin(1, 15, 1, 1))

pdf(paste0(fig_path, "/or_reri.pdf"), width = 8, height = 3)
or_plot + reri_plot #+ plot_layout(widths = c(3,3))
dev.off()

# also get these values in a table
n_digits <- 3
res_tab <- res_df %>%
  mutate(model = str_remove(mod.name, ".RData"),
         or_ci = paste0(formatC(round(point.est, n_digits), format = "f", digits = n_digits),
                        " (", formatC(round(`2.5%`, n_digits), format = "f", digits = n_digits),
                        ", ", formatC(round(`97.5%`, n_digits), format = "f", digits = n_digits),
                        ")"),
         pval = ifelse(pval < 0.001, "<0.001", 
                       formatC(round(pval, n_digits), format = "f", digits = n_digits)))

res_tab$estimand <- rownames(res_df)

res_tab <- res_tab %>%
  select(estimand, or_ci, pval) %>%
  mutate(estimand = factor(estimand, levels = c("OR10", "OR01", "OR11", "RERI"))) %>%
  arrange(estimand)

# save this table (make a table with all model results in 5_OR_table.R)
saveRDS(res_tab, file = paste0("data/intermediate/OR_tables/OR_table_", dataset, 
        "_", str_remove(mod.name, ".RData"), ".rds"))


#---------- visualize effect of just temp ----------#

# how does OR10 change for different temps when we hold PM at the median?
# OR 10 calc doesn't require bad PM value --> just good, which will be median

# reference temp
ref_temp <- med_temp
#ref_temp <- 60

# vector of comparison temps
if(dataset %in% c("tmax_temp", "hi_temp")){
  comp_temp_vec <- seq(60, 115, by = 1)
} else {
  comp_temp_vec <- seq(0.3, 1, by = 0.001)
}

# constant PM2.5
constant_pm25 <- med_pm25

# initialize data frame for results
temp_res <- as.data.frame(matrix(nrow = length(comp_temp_vec), ncol = 4)) |>
  setNames(c("comp_expo", "OR", "lower", "upper"))

temp_res$comp_expo <- comp_temp_vec


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
  ci95 <- quantile(OR, c(0.05, 0.95))
  
  temp_res[j,2:4] <- c(point.est, ci95)
  
}


#---------- visualize effect of just PM2.5 ----------#

# reference PM2.5
ref_pm25 <- med_pm25
#ref_pm25 <- 5

# vector of comparison PM2.5
if(str_detect(mod.name, "lowPM")){
  comp_pm25_vec <- seq(0, 18, by = 1)
}else{
  comp_pm25_vec <- seq(0, 100, by = 1)
}


# constant temp
constant_temp <- med_temp

# initialize data frame for results
pm25_res <- as.data.frame(matrix(nrow = length(comp_pm25_vec), ncol = 4)) |>
  setNames(c("comp_expo", "OR", "lower", "upper"))

pm25_res$comp_expo <- comp_pm25_vec


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


# temp plot
temp_effect_plot <- temp_res |>
  ggplot(aes(x = comp_expo, y = OR)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_line(col = "#ff9e00") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#ff9e00") +  
  #geom_vline(xintercept = med_temp, lty = 2, col = "gray40") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(#breaks = c(0.8, 0.9, 1, 1.1, 1.2),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.77, 1.23)) +
  labs(#title = mod.name,
       x = "Daily maximum temperature (\u00B0F)",
       y = "OR (compared to median)") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.2)
  )

# pm25 plot
pm25_effect_plot <- pm25_res |>
  ggplot(aes(x = comp_expo, y = OR)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_line(col = "purple3") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "purple3") +
  #geom_vline(xintercept = med_pm25, lty = 2, col = "gray40") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(#breaks = c(0.95, 1, 1.05),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.92, 1.08)) +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       #y = expression("OR (compared to median " * PM[2.5] * ")")
       y = "OR (compared to median)") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.2)
        )

# align plots (so the plotting space is the same)
temp_pm25_plots <- cowplot::align_plots(temp_effect_plot, pm25_effect_plot, align = "hv")

# plot
pdf(paste0(fig_path, "independent_effects.pdf"), height = 3, width = 6.5)
plot_grid(temp_pm25_plots[[1]], temp_pm25_plots[[2]])
dev.off()

#---- instead, try using facet_grid
# I don't like this because they need separate X labels

# temp_res$exposure <- "temp"
# pm25_res$exposure <- "pm25"
# 
# temp_res$expo_median <- med_temp
# pm25_res$expo_median <- med_pm25
# 
# temp_pm25_res <- rbind(temp_res, pm25_res)
# temp_pm25_res$exposure <- factor(temp_pm25_res$exposure, levels = c("temp", "pm25"))
# 
# # plot colors
# expo_cols <- c("temp" = "purple3", "pm25" = "orange2")
# 
# # facet labels
# expo_labs <- c("temp" = "Daily~max~temperature~(degree*F)", 
#                "pm25" = "Three-day~PM[2.5]~(mu*g/m^3)")
# 
# # plot both
# temp_pm25_res |>
#   ggplot(aes(x = comp_expo, y = OR, col = exposure, fill = exposure)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA) +
#   geom_hline(yintercept = 1) +
#   geom_vline(aes(xintercept = expo_median), lty = 2, col = "gray40") +
#   labs(x = "Comparison exposure level",
#        y = "OR (compared to median)") +
#   scale_color_manual(values = expo_cols) +
#   scale_fill_manual(values = expo_cols) +
#   facet_wrap(~exposure, 
#              scales = c("free"), 
#              labeller = as_labeller(expo_labs, default = label_parsed)) +
#   theme_minimal() +
#   theme(legend.position = "none",
#         panel.spacing = unit(2, "lines"))


#---------- risk surface plot ----------#

# # this only works is using a tensor (slow)
# plot(conditional_smooths(stan.mod), ask = FALSE)

# try creating plot manually

# get pairs of PM and temperature values 
temp_vals <- seq(60, 115, length = 10)
if(str_detect(mod.name, "lowPM")){
  pm25_vals <- seq(0, 18, length = 10)
}else{
  pm25_vals <- seq(0, 100, length = 10)
}
temp_pm25_grid <- expand.grid(temp = temp_vals, pm25 = pm25_vals)
temp_pm25_grid$OR11 <- NA

# loop through temperature and PM2.5 pairs
for(i in 1:nrow(temp_pm25_grid)){
  
  # get current temp
  temp_i <- temp_pm25_grid$temp[i]
  pm25_i <- temp_pm25_grid$pm25[i]
  
  # get each basis for current values
  cf11 <- get_cf_basis(temp = temp_i, pm25 = pm25_i)
  cf00 <- get_cf_basis(temp = med_temp, pm25 = med_pm25)
  
  # initialize vector for odds ratios
  # comparing exposure to both current values to exposure to both at median
  OR11 <- c()
  
  # loop through beta.post to get OR11 estimates
  for(j in 1:nrow(beta.post)){
    # get odds ratio comparing 11 to reference
    OR11[j] <- exp(cf11 %*% beta.post[j,] - cf00 %*% beta.post[j,])
  }
  
  # get posterior mean
  temp_pm25_grid$OR11[i] <- mean(OR11)
    
}


# plot risk surface
fill_scale <- scale_fill_gradient2(low = "#e4cd05", mid = "white", high = "#b11226",
                                   midpoint = 1,
                                   limits = c(0.8, 1.3)
                                   )
pdf(paste0(fig_path, "risk_surface.pdf"), width = 4, height = 4.5)
temp_pm25_grid |>
  ggplot() + 
  geom_tile(aes(x = temp, y = pm25, fill = OR11)) +
  labs(x = "Daily maximum temperature (\u00B0F)",
       y = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       #fill = expression("OR (vs. median temperature & " * PM[2.5] * ")")
       #fill = "Odds ratio (compared to median exposures)"
       fill = "Odds ratio (compared to \u00d7)"
       ) +
  fill_scale +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(x = med_temp, y = med_pm25, shape = 4) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 15,
                               barheight = 0.5,
                               ticks.colour = NA))
dev.off()


