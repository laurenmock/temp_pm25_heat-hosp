###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Make plots/tables with results from multiple models
# Input: 
#   Processed results from script 4
# Outputs: 
#   Figures and tables with results from multiple models
###############################################################

library(xtable)
library(tidyverse)
library(ggplot2)
library(patchwork)

# path to tables for each model
path_or <- "data/intermediate/OR_tables/"

# path for saving figures
fig_path <- paste0("results/figures/")

# load tables from each model
res_main <- readRDS(file = paste0(path_or, "OR_table_tmax_temp_spline_trim95.rds"))
res_temp3day <- readRDS(file = paste0(path_or, "OR_table_tmax_temp_spline_trim95_temp3day.rds"))
res_trim99 <- readRDS(file = paste0(path_or, "OR_table_tmax_temp_spline_trim99.rds"))
res_t2 <- readRDS(file = paste0(path_or, "OR_table_tmax_temp_t2_trim95.rds"))

# join results from all models
res_all <- rbind(res_main, res_temp3day, res_trim99, res_t2)


#--------------------------- COMPARISON TABLE ---------------------------#

n_digits_or <- 2
res_tab <- res_all %>%
  mutate(or_ci = paste0(formatC(round(point.est, n_digits_or), format = "f", digits = n_digits_or),
                        " (", formatC(round(`2.5%`, n_digits_or), format = "f", digits = n_digits_or),
                        ", ", formatC(round(`97.5%`, n_digits_or), format = "f", digits = n_digits_or),
                        ")"))

res_tab <- res_tab %>%
  select(model, estimand, or_ci)

res_tab <- pivot_wider(res_tab,
            id_cols = c("estimand"),
            names_from = "model",
            values_from = c("or_ci"))

print(xtable(res_tab,
             type = latex,
             label = "tab:or_results",
             caption = ""),
      file = "results/tables/or_reri_results.tex",
      sanitize.text.function = identity,
      include.rownames = FALSE)


#--------------------------- COMPARISON FIGURE ---------------------------#

# set levels
res_all$estimand <- factor(res_all$estimand,
                          levels = rev(c("RERI",
                                         "OR11",
                                         "OR01",
                                         "OR10")))
# column for OR/RERI facet
res_all <- res_all %>%
  mutate(group = ifelse(estimand == "RERI", "Interaction", "Odds ratios"),
         null_val = ifelse(estimand == "RERI", 0, 1))
res_all$group <- factor(res_all$group, 
                         levels = c("Odds ratios", "Interaction"))

# plot OR
or_plot <- res_all %>%
  filter(estimand != "RERI") %>%
  ggplot(aes(x = estimand, y = point.est, group = model, col = model)) + 
  geom_hline(yintercept = 1, col = "gray60", lty = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                width = 0, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "Odds ratio", color = "Model") +
  scale_color_manual(values = c("spline_trim95" = "#26547C",
                                "spline_trim95_temp3day" = "brown3",
                                "spline_trim99" = "goldenrod",
                                "t2_trim95" = "dodgerblue"),
                     labels = c("spline_trim95" = "Main",
                                "spline_trim95_temp3day" = "Three-day temperature exposure",
                                "spline_trim99" = expression(PM[2.5]*" trimmed at 99th percentile"),
                                "t2_trim95" = "Tensor product")) +
  scale_x_discrete(labels = c("OR10" = "Temperature",
                              "OR01" = expression("PM"["2.5"]),
                              "OR11" = expression("Temperature and PM"["2.5"]))) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.background = element_blank())

# plot RERI
reri_plot <- res_all %>%
  filter(estimand == "RERI") %>%
  ggplot(aes(x = estimand, y = point.est, group = model, col = model)) + 
  geom_hline(yintercept = 0, col = "gray60", lty = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                width = 0, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "RERI", color = "Model") +
  scale_color_manual(values = c("spline_trim95" = "#26547C",
                                "spline_trim95_temp3day" = "brown3",
                                "spline_trim99" = "goldenrod",
                                "t2_trim95" = "dodgerblue"),
                     labels = c("spline_trim95" = "Main",
                                "spline_trim95_temp3day" = "Three-day temperature exposure",
                                "spline_trim99" = expression(PM[2.5]*" trimmed at 99th percentile"),
                                "t2_trim95" = "Tensor product")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# align OR and RERI plots
(or_plot | reri_plot) + 
  plot_layout(widths = c(2, 1),
              guides = "collect") & 
  theme(legend.position = "bottom",  
        legend.box = "vertical")

# save plot
ggsave(paste0("results/figures/OR_RERI_across_models.pdf"), height = 4, width = 8)



#--------------------------- INDEPENDENT EFFECTS ---------------------------#

# results from the three sensitivity analyses

path_indep <- "data/intermediate/independent_effects/"

# load plotting data from script 4
temp_temp3day <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_temp_tmax_temp_spline_trim95_temp3day.rds"))
pm25_temp3day <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_pm25_tmax_temp_spline_trim95_temp3day.rds"))
temp_trim99 <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_temp_tmax_temp_spline_trim99.rds"))
pm25_trim99 <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_pm25_tmax_temp_spline_trim99.rds"))
temp_t2 <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_temp_tmax_temp_t2_trim95.rds"))
pm25_t2 <- readRDS(file = paste0(path_indep, 
                                       "independent_effects_pm25_tmax_temp_t2_trim95.rds"))

# bind temperature results together across models
indep_temp <- rbind(temp_temp3day, temp_trim99, temp_t2)

# bind PM2.5 results together across models
indep_pm25 <- rbind(pm25_temp3day, pm25_trim99, pm25_t2)

# get facet labels
model_labs <- c(
  spline_trim95_temp3day = '"Three-day temperature exposure"',
  spline_trim99 = 'PM[2.5]*" trimmed at 99th percentile"',
  t2_trim95 = '"Tensor product"'
)

# temp plot
temp_plot <- indep_temp |>
  ggplot(aes(x = comp_expo, y = OR)) +
  geom_hline(yintercept = 1, col = "black", lty = 1, linewidth = 0.3) +
  geom_line(col = "#ffad4d") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#ffad4d") +  
  geom_vline(xintercept = 29.6, lty = 2, col = "gray40") + # median
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.82, 1.25)) +
  labs(x = "Daily maximum temperature (\u00B0C)",
       y = "Odds ratio (relative to median)") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        strip.text = element_text(size = 12)) +
  facet_grid(~model,
             labeller = labeller(model = model_labs, .default = label_parsed))

# pm25 plot
pm25_plot <- indep_pm25 |>
  ggplot(aes(x = comp_expo, y = OR)) +
  geom_hline(yintercept = 1, col = "black", lty = 1, linewidth = 0.3) +
  geom_line(col = "#631ea2") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#631ea2") +
  geom_vline(xintercept = 8.9, lty = 2, col = "gray40") + # median
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.84, 1.25)) +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       y = "Odds ratio (relative to median)") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        strip.text = element_text(size = 12)) +
  facet_grid(~model,
             labeller = labeller(model = model_labs, .default = label_parsed))

# align plots
temp_plot / pm25_plot

# save
ggsave(paste0(fig_path, "sensitivity_independent_effects.pdf"), width = 9, height = 6)


#--------------------------- RISK SURFACE ---------------------------#

# results from the three sensitivity analyses

path_surface <- "data/intermediate/risk_surface/"

# load plotting data from script 4
surface_temp3day <- readRDS(file = paste0(path_surface, "risk_surface_tmax_temp_spline_trim95_temp3day.rds"))
surface_trim99 <- readRDS(file = paste0(path_surface, "risk_surface_tmax_temp_spline_trim99.rds"))
surface_t2 <- readRDS(file = paste0(path_surface, "risk_surface_tmax_temp_t2_trim95.rds"))

# bind results together
surface <- rbind(surface_temp3day, surface_trim99, surface_t2)

# color scale
fill_scale <- scale_fill_gradient2(low = "#056608", mid = "white", high = "#a73672",
                                   midpoint = 1,
                                   limits = c(0.84, 1.331)
)

# plot
surface |>
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
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 15,
                               barheight = 0.5,
                               ticks.colour = NA)) + 
  facet_grid(~model,
             labeller = labeller(model = model_labs, .default = label_parsed))

# save
ggsave(paste0(fig_path, "sensitivity_risk_surface.pdf"), width = 8.5, height = 4)
