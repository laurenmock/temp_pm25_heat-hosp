###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Understand data
# Input: 
#   Case-crossover data
# Output: 
#   Summary statistics/EDA
###############################################################

library(data.table)
library(ggplot2)
library(table1)
library(ggpubr)
library(scales)
library(tidyverse)
library(sf)
library(xtable)
library(haven)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(ggpointdensity)


#---------- load case-crossover data ----------#

# Load trimmed case-crossover data
load(paste0("data/intermediate/cco_trim95_tmax_temp.RData"))
cco <- cco_trim95
rm(cco_trim95); gc()


#---------- convert temperature to Celsius ----------#

cco[, temp_lag0 := (temp_lag0 - 32) * (5/9)]


#---------- load zip to county data ----------#

zip2county <- read.csv("data/raw/zip2county_master_xwalk_2010_2023_tot_ratio_one2one.csv")

# select columns of interest
zip2county <- zip2county %>%
  filter(year == "2010") %>%
  rename("zipcode" = "zip") %>%
  select(zipcode, county)

# make cco dataset zipcode integer
cco <- cco %>%
  mutate(zipcode = as.integer(zipcode))

# merge with dataset
cco <- merge (cco, zip2county, by = "zipcode", all.x = TRUE, all.y = FALSE)


#---------- Basics ----------#

# Number of individuals/first hospitalizations
length(unique(cco[,bene_id]))

# Age distribution
cco[cases == 1] |>
  ggplot() +
  geom_histogram(aes(age), 
                 #col = "white", 
                 fill = "#a0c293", 
                 #alpha = 0.5,
                 binwidth = 1) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Age on case date", y = "Number of individuals") +
  theme_minimal()

# Number of control days after trimming
cco[cases == 1] |>
  ggplot() +
  geom_bar(aes(as.factor(num_control_trim)), 
           fill = "#a0c293") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Number of control days", y = "Number of individuals") +
  theme_minimal()

# Case year
cco[cases == 1] |>
  ggplot() +
  geom_bar(aes(caseyr), col = "white", fill = "#a0c293") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Case year", y = "Number of Individuals") +
  theme_minimal()

# Case dow
cco[cases == 1] |>
  ggplot() +
  geom_bar(aes(dow), col = "white", fill = "#a0c293") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Case day of week", y = "Number of Individuals") +
  theme_minimal()

# Case month
cco <- cco |> 
  mutate(month_text = case_when(
    month == 6 ~ "June",
    month == 7 ~ "July",
    month == 8 ~ "August",
    month == 9 ~ "September"
  ))
cco$month_text <- factor(cco$month_text,
                         levels = c("June", "July", "August", "September"))
cco[cases == 1] |>
  ggplot() +
  geom_bar(aes(month_text), col = "white", fill = "#a0c293") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Case month", y = "Number of Individuals") +
  theme_minimal()

# Control month
cco[cases == 0] |>
  ggplot() +
  geom_bar(aes(month_text), col = "white", fill = "#a0c293") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Control month", y = "Number of Control dates") +
  theme_minimal()


#---------- heat and PM2.5 exposures ----------#

# temperature (case days)
temp_dist <- cco |>
  filter(cases == 1) |>
  ggplot() +
  geom_histogram(aes(temp_lag0),
                 fill = "#f5ca7b", 
                 #alpha = 0.5,
                 binwidth = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = label_comma(), 
                     expand = c(0, 0),
                     limits = c(0, 20000)) +
  geom_vline(xintercept = median(cco$temp_lag0), lty = 2, color = "gray30") +
  geom_vline(xintercept = quantile(cco$temp_lag0, 0.95), lty = 2, color = "gray30") +
  labs(x = "Daily maximum temperature (\u00B0C)", 
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.2)
  )

# PM2.5 over three days (case days)
pm25_dist <- cco |>
  filter(cases == 1) |>
  ggplot() +
  geom_histogram(aes(pm25_3day),
                 fill = "#c0aada", 
                 binwidth = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = label_comma(), 
                     expand = c(0, 0),
                     limits = c(0, 20000)) +
  geom_vline(xintercept = median(cco$pm25_3day), lty = 2, color = "gray30") +
  geom_vline(xintercept = quantile(cco$pm25_3day, 0.95), lty = 2, color = "gray30") +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"), 
       y = "") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.2)
  )

# align plots (so the plotting space is the same)
temp_pm25_plots <- cowplot::align_plots(temp_dist, pm25_dist, align = "hv")

# plot
plot_grid(temp_pm25_plots[[1]], temp_pm25_plots[[2]])
ggsave("results/figures/temp_pm25_dist.pdf", width = 7, height = 3)

# get range of temp and pm2.5 values on case date
cases <- cco %>%
  filter(cases == 1)

range(cases$temp_lag0)
quantile(cases$temp_lag0, c(0.05, 0.95))
median(cases$temp_lag0)
mean(cases$temp_lag0)

range(cases$pm25_3day) # after trimming (see bottom of this script for max before trimming)
median(cases$pm25_3day)


#---------- boxplots of exposure on case vs control days ----------#

cco %>%
  mutate(cases = ifelse(cases == 1, "Case", "Control")) %>%
  ggplot() +
  geom_boxplot(aes(temp_lag0, x = as.factor(cases)), fill = "#f5ca7b") +
  labs(x = "", y = "Daily max. temperature (\u00B0C)") +
  theme_minimal()

cco %>%
  mutate(cases = ifelse(cases == 1, "Case", "Control")) %>%
  ggplot() +
  geom_boxplot(aes(pm25_3day, x = as.factor(cases)), fill = "#c0aada") +
  labs(x = "", y = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")")) +
  theme_minimal()


#---------- temp vs. PM2.5 ----------#

cases |>
  ggplot(aes(x = temp_lag0, y = pm25_3day)) +
  geom_bin2d() + # from ggpointdensity package
  labs(x = "Daily maximum temperature (\u00B0C)",
       y = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       fill = "Number of case days") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis(option = "A") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 20,
                               barheight = 0.5,
                               ticks.colour = NA))
ggsave("results/figures/temp_vs_pm25.pdf", height = 5.5, width = 5)


#---------- maps (exposures and outcome) ----------#

# set up data for plotting (while cco is a data table--faster)
plot_df <- cco %>%
  filter(cases == 1) %>%
  group_by(county) %>%
  summarise(mean_temp_lag0 = mean(temp_lag0),
            mean_pm25_3day = mean(pm25_3day),
            n = n())

# load at risk people in our Medicare data by county
at_risk <- read_sas(paste0("data/raw/atrisk_cty.sas7bdat"))

# join with plot_df
plot_df <- left_join(plot_df, at_risk, by = "county")

# new column of interest: hosp rate
plot_df <- plot_df %>%
  mutate(hosp_rate = n/at_risk * 1000)


#---- get county/state geometry

# load state shapefile
state_sf <- read_sf("data/raw/shapefiles_state/cb_2018_us_state_20m.shp") %>%
  filter(!(STATEFP %in% c("02", "15", "66", "72", "60", "69", "78")))

# load county shapefile
county_sf <- read_sf("data/raw/shapefiles_county/cb_2018_us_county_20m.shp") %>%
  filter(!(STATEFP %in% c("02", "15", "66", "72", "60", "69", "78"))) %>%
  mutate(county = as.numeric(GEOID))

# join county shapefile with plotting data
plot_df <- full_join(county_sf, plot_df, by = c("county"), )


##### temp map #####
temp_map <- plot_df %>%
  mutate(mean_temp_lag0 = ifelse(n > 10, 
                                 mean_temp_lag0, NA)) %>% # filter to > 10 cases
  mutate(temp_cat = cut(mean_temp_lag0, 
                        right = FALSE,
                        breaks = c(17, 27, 29, 32, 41),
                        labels = c("17\u201327", "27\u201329",
                                   "29\u201332", "32\u201341")
                        )) %>%
  ggplot() +
  geom_sf(aes(fill = temp_cat), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "black") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = "Daily maximum temperature (\u00B0C)") +
  scale_fill_manual(values = viridis(4, option = "D"),
                    na.translate = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.ticks = element_blank())


##### pm25 map #####
pm25_map <- plot_df %>%
  mutate(mean_pm25_3day = ifelse(n > 10, 
                                 mean_pm25_3day, NA)) %>% # filter to > 10 cases
  mutate(pm25_cat = cut(mean_pm25_3day, 
                        right = FALSE,
                        breaks = c(2, 7, 9, 11, 15),
                        labels = c("3\u20137", "7\u20139",
                                   "9\u201311", "11\u201314")
  )) %>%
  ggplot() +
  geom_sf(aes(fill = pm25_cat), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "black") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")")) +
  scale_fill_manual(values = viridis(4, option = "D"),
                    na.translate = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.ticks = element_blank())


###### rate of hospitalization map #####

hosp_map <- plot_df %>%
  mutate(hosp_rate = ifelse(n > 10,
                            hosp_rate, NA)) %>% # filter to > 10 cases
  mutate(hosp_cat = cut(hosp_rate,
                        right = FALSE,
                        breaks = c(1, 4, 7, 10, 59),
                        labels = c("1\u20134", "4\u20137",
                        "7\u201310", "10\u201358")
  )) %>%
  ggplot() +
  geom_sf(aes(fill = hosp_cat), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "black") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = "Heat-related hospitalizations \nper 1,000 people at risk") +
  scale_fill_manual(values = viridis(4, option = "D"),
                    na.translate = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.ticks = element_blank())


# align plots (so the plotting space is the same)
all_maps <- cowplot::align_plots(temp_map, pm25_map, hosp_map, align = "hv")

# plot
plot_grid(all_maps[[1]], all_maps[[2]], all_maps[[3]], nrow = 1)
ggsave("results/figures/temp_pm25_hosp_maps.pdf", height = 3, width = 10)


##############################################################

# dataset before trimming 
load(paste0("data/intermediate/cco_full_tmax_temp.RData"))

# get 95th and 99th percentiles of 3-day PM2.5
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95) # 18.4 ug/m3
pm25_3day_q99 <- quantile(cco$pm25_3day, 0.99) # 23.7 ug/m3

# get dist of PM2.5 before trimming
cco |>
  filter(cases == 1) |>
  ggplot() +
  geom_histogram(aes(pm25_3day),
                 fill = "#c0aada", 
                 binwidth = 0.5) +
  geom_vline(xintercept = pm25_3day_q95, col = "black") +
  geom_vline(xintercept = pm25_3day_q99, col = "black") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150)) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       y = "Frequency") +
  theme_minimal()
ggsave("results/figures/pm25_dist_pretrim.pdf", width = 8, height = 3)

# get range of PM2.5 values (before trimming) on case dates
cases <- cco %>%
  filter(cases == 1)

max(cases$pm25_3day)
max(cco$pm25_3day)

