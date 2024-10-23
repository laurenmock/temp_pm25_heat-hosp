###############################################################
# Heat-related outcomes
# Author: Lauren Mock
# Aim: Understand data
# Input: Case-crossover data
# Output: Summary statistics/EDA
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

#---------- user inputs ----------#

# choose temperature metric
dataset <- "tmax_temp"
#dataset <- "tmax_pcnt"
#dataset <- "hi_temp"
#dataset <- "hi_pcnt"


#---------- load case-crossover data ----------#

# Load case-crossover data
load(paste0("data/intermediate/cco_trimmed_", dataset, ".RData"))

# open PDF to save plots
#pdf("results/EDA.pdf")


#---------- load zip to county data ----------#

zip2county <- read.csv("data/raw/zip2county_master_xwalk_2010_2023_tot_ratio_one2one.csv")

# select columns of interest
zip2county <- zip2county %>%
  filter(year == "2010") %>%
  rename("zipcode" = "zip") %>%
  select(zipcode, county)

# make cco dataset zipcode integer (need to check this is working properly)
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


#--- Med use
# what % are not taking any meds?
# sum(is.na(cco[cases == 1, user_Antichol])) / nrow(cco[cases == 1])

# # what % are taking each of the three meds?
# cco |>
#   # filter to people taking at least one med
#   filter(cases == 1) |>
#   select(user_Antichol, user_Stim, user_LoopD) |>
#   colMeans(na.rm = TRUE)
#   #gather(key = "key", value = "value") |>
#   # group_by(key) |>
#   # summarise(frequency = n_distinct(value))

#test[,table(value), by = key]


# cco[cases == 1] |>
#   group_by(steroid_use) |>
#   summarise(count = n()) |>
#   mutate(percent = 100 * count/sum(count)) |>
#   ggplot() +
#   geom_bar(aes(x = steroid_use, y = percent), 
#            col = "white", fill = "#a0c293",
#            stat = "identity") +
#   labs(x = "Steroid use on case and control days", y = "% of individuals") +
#   scale_y_continuous(labels = label_comma()) +
#   theme_minimal()

# # steroid use on control days before, steroid use on case day, steroid use on control days after
# cco |>
#   group_by(control_time, steroid_lag0) |>
#   summarise(count = n()) |>
#   mutate(percent = 100 * count/sum(count)) |>
#   ggplot() +
#   geom_bar(aes(x = control_time, y = percent, fill = factor(steroid_lag0)), 
#            col = "white", stat = "identity") +
#   labs(x = "", y = "% of days", fill = "") +
#   scale_y_continuous(labels = label_comma()) +
#   scale_fill_manual(values = c(`0` = "skyblue", 
#                                    `1` = "darkorange"),
#                         labels = c("Not taking steroids", "Taking steroids")) +
#   theme_minimal()



#---------- look at inhosp variable ----------#

# look at inhosp variable
100 * nrow(cco[cases == 0 & inhosp == 1])/nrow(cco[cases == 0]) # 8%
# leaving these for now, but we could remove them

# plot.new()
# text(x = 0.5, y = 0.95, paste0(round(pct_inhosp, 1), 
#                                "% of control days were during hospital stays"))


#---------- heat and PM2.5 exposure ----------#

# cco |>
#   ggplot() +
#   geom_density(aes(temp_lag0, col = as.factor(cases))) +
#   theme_minimal()
# 
# cco |>
#   ggplot() +
#   geom_density(aes(pm25_lag0, col = as.factor(cases))) +
#   theme_minimal()

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
                     limits = c(0, 10500)) +
  geom_vline(xintercept = median(cco$temp_lag0), lty = 2, color = "gray30") +
  geom_vline(xintercept = quantile(cco$temp_lag0, 0.95), lty = 2, color = "gray30") +
  labs(x = "Daily maximum temperature (\u00B0F)", 
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
                 binwidth = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = label_comma(), 
                     expand = c(0, 0),
                     limits = c(0, 10500)) +
  geom_vline(xintercept = median(cco$pm25_3day), lty = 2, color = "gray30") +
  geom_vline(xintercept = quantile(cco$pm25_3day, 0.95), lty = 2, color = "gray30") +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"), 
       y = "") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.2)
  )

# align plots (so the plotting space is the same)
temp_pm25_plots <- align_plots(temp_dist, pm25_dist, align = "hv")

# plot
pdf("results/figures/temp_pm25_dist.pdf", width = 7, height = 3)
plot_grid(temp_pm25_plots[[1]], temp_pm25_plots[[2]])
dev.off()

# get range of temp and pm values on case date
cases <- cco %>%
  filter(cases == 1)

range(cases$temp_lag0)
range(cases$pm25_3day) # after trimming (see bottom of this script for max before trimming)
median(cases$pm25_3day)

#---------- boxplots of exposure on case vs control days ----------#

cco %>%
  mutate(cases = ifelse(cases == 1, "Case", "Control")) %>%
  ggplot() +
  geom_boxplot(aes(temp_lag0, x = as.factor(cases)), fill = "#f5ca7b") +
  labs(x = "", y = "Daily max. temperature (\u00B0F)") +
  theme_minimal()

cco %>%
  mutate(cases = ifelse(cases == 1, "Case", "Control")) %>%
  ggplot() +
  geom_boxplot(aes(pm25_3day, x = as.factor(cases)), fill = "#c0aada") +
  labs(x = "", y = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")")) +
  theme_minimal()


##################################################

# # different temp exposures
# 
# dataset <- "tmax_temp"
# #dataset <- "tmax_pcnt"
# #dataset <- "hi_temp"
# #dataset <- "hi_pcnt"
# 
# load(paste0("data/intermediate/cco_", dataset, ".RData"))
# 
# # Age distribution
# cco[cases == 1] |>
#   ggplot() +
#   geom_histogram(aes(temp_lag0), col = "white", fill = "#a0c293",
#                  binwidth = 3
#                 ) +
#   scale_y_continuous(labels = label_comma()) +
#   labs(x = paste0(dataset, " on case date"), y = "Count") +
#   theme_minimal()
# 


##################################################
##############    MAPS  ##########################
##################################################

# # load zip to state file
# zip_to_state <- read_sas(paste0("data/raw/zipcode_to_state.sas7bdat"))
# 
# # merge with patient data
# cco <- left_join(cco, zip_to_state, by = c("zipcode"))

#--- get all the plotting data I want while cco is a data table (faster)

# # get data to plot--mean case day temp in each state
# plot_df <- cco %>%
#   filter(cases == 1) %>%
#   group_by(State) %>%
#   summarise(mean_temp_lag0 = mean(temp_lag0),
#             mean_pm25_3day = mean(pm25_3day),
#             n = n())

plot_df <- cco %>%
  filter(cases == 1) %>%
  group_by(county) %>%
  summarise(mean_temp_lag0 = mean(temp_lag0),
            mean_pm25_3day = mean(pm25_3day),
            n = n())


#--- now merge and plot

# load state shapefile
state_sf <- read_sf("data/raw/shapefiles_state/cb_2018_us_state_20m.shp") %>%
  filter(!(STATEFP %in% c("02", "15", "66", "72", "60", "69", "78")))

# load county shapefile
county_sf <- read_sf("data/raw/shapefiles_county/cb_2018_us_county_20m.shp") %>%
  filter(!(STATEFP %in% c("02", "15", "66", "72", "60", "69", "78"))) %>%
  mutate(county = as.numeric(GEOID))

# join shapefile with plotting data
#plot_df <- right_join(county_sf, plot_df, by = c("STUSPS" = "State"))
plot_df <- full_join(county_sf, plot_df, by = c("county"), )

# load at risk people in Medicare
at_risk <- read_sas(paste0("data/raw/atrisk_cty.sas7bdat"))

# join with plot_df
plot_df <- left_join(plot_df, at_risk, by = "county")




#color_pal_temp <- c("#ffffd2", "#FFFFB2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")
color_pal_temp <- c("#ffffd2", "#FFFFB2", "#fecc5c", "#fd8d3c", "#cc5500", "#7b3804")

# temp map
temp_map <- plot_df %>%
  mutate(mean_temp_lag0 = ifelse(n > 10, mean_temp_lag0, NA)) %>% # filter to > 10 cases
  ggplot() +
  geom_sf(aes(fill = mean_temp_lag0), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "gray90") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = "Mean daily max. temperature (\u00B0F)") +
  #scale_fill_gradient(low = "#fff7ed", high = "#ff6600") +
  # scale_fill_gradientn(name = waiver(),
  #                      low = "#ffffbf", 
  #                      high = "red3") +
  #scale_fill_gradientn(colors = brewer.pal(n = 8, name = "Oranges")) +
  #scale_fill_viridis(option = "F", na.value = "gray70", direction = -1) +
  scale_fill_gradientn(colors = color_pal_temp, na.value = "gray70") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 10,
                               barheight = 0.5))

color_pal_pm25 <- c("#ede7f6", "#b19cd9", "#7e5bca", "#3f007d", "#130a2e")

# pm25 map
pm25_map <- plot_df %>%
  mutate(mean_pm25_3day = ifelse(n > 10, mean_pm25_3day, NA)) %>% # filter to > 10 cases
  ggplot() +
  geom_sf(aes(fill = mean_pm25_3day), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "gray90") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = expression("Mean three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")")) +
  #scale_fill_gradient(low = "#e9e6f3", high = "#692d94", na.value = "gray90") +
  #scale_fill_viridis(option = "D", na.value = "gray70") +
  scale_fill_gradientn(colors = color_pal_pm25, na.value = "gray70") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 10,
                               barheight = 0.5))

# rate of hosp map 

# # further aggregate plot data by state
# plot_state_df <- plot_df %>%
#   group_by(STATEFP) %>%
#   summarise(n = sum(n)) %>%
#   st_drop_geometry()

# # join shapefile with plotting data
# plot_state_sf <- full_join(state_sf, plot_state_df, by = c("STATEFP"))
# 
# # merge with plot data
# plot_state_df <- left_join(plot_state_df, at_risk, by = c("STATEFP" = "State"))

# new column of interest
plot_df <- plot_df %>%
  mutate(hosp_rate = n/at_risk * 1000)

color_pal_atrisk <- colorRampPalette(c(#"#f9eaed", 
                                       "#f2dede", "#ecc0c9", "#e079a6", "#a2004c", "#750d3e", "#531f20"))

atrisk_map <- plot_df %>%
  mutate(hosp_rate = ifelse(n > 10, hosp_rate, NA)) %>% # filter to > 10 cases
  mutate(hosp_rate = ifelse(hosp_rate > 40, 40, hosp_rate)) %>% # cap max val
  ggplot() +
  geom_sf(aes(fill = hosp_rate), col = NA) +
  geom_sf(data = state_sf, fill = NA, col = "gray90") +
  coord_sf(crs = st_crs(5070)) +
  labs(fill = "Heat-related hosp. / 1,000 people at risk") +
  scale_fill_gradientn(colors = color_pal_atrisk(10), na.value = "gray70",
                       breaks = c(10, 20, 30, 40),
                       labels = c("10", "20", "30", "40+")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 10,
                               barheight = 0.5))

pdf("results/figures/temp_pm25_hosp_maps.pdf", height = 3, width = 9)
ggarrange(temp_map, pm25_map, atrisk_map, nrow = 1)
dev.off()


# need to figure out which counties to remove --> don't worry about it
# # get total # at risk
# at_risk <- at_risk %>%
#   filter(!State %in% c("AA", "AR", "AK", "AS", "AE", "AP", "GU", "HI", "MP", "PR", "DC"))
# total_at_risk <- sum(at_risk$at_risk)


# #--- number of cases by zip code (or maybe by state)
# 
# 
# #--- some kind of maps of PM2.5 and temperature 
# 
# # I don't even have this for the whole study period (only case/controls and lags)
# 
# 
# # load state and zcta shapefiles
# state_sf <- read_sf("data/raw/shapefiles_state/cb_2018_us_state_20m.shp") %>%
#   filter(!(STATEFP %in% c("02", "15", "66", "72", "60", "69", "78")))
# zcta_sf <- read_sf("data/raw/shapefiles_zcta/cb_2018_us_zcta510_500k.shp")
# 
# 
# # # FILTER STATES MORE FOR NOW
# # state_sf <- state_sf %>%
# #   filter(STATEFP %in% c("05", "29"))
# 
# 
# # # try mapping mean temp_lag0 on case days (SLOW)
# # map_temp_df <- cco_sf %>%
# #   group_by(ZCTA5CE10) %>%
# #   summarise(mean_temp_lag0 = first(temp_lag0))
# 
# 
# # and maybe I should do all of this for counties instead of states? but borders get messy
# # could also just do zip codes, but that doesn't include a lot of areas
# 
# # get zip/state intersections (slow--use two states only for now)
# zip_state_intersect <- st_intersection(zcta_sf, state_sf)
# # warning is fine
# 
# # new column to get intersection area
# zip_state_intersect <- zip_state_intersect %>%
#   mutate(intersection_area = st_area(geometry))
# 
# # # check 65733
# # zip_state_intersect %>%
# #   filter(ZCTA5CE10 == "65733") %>%
# #   pull(intersection_area)
# 
# # for each county, filter to state with largest area
# zip_to_state <- zip_state_intersect %>%
#   group_by(ZCTA5CE10) %>%
#   filter(intersection_area == max(intersection_area))
# 
# # # check 65733 again
# # zip_to_state %>%
# #   filter(ZCTA5CE10 == "65733") %>%
# #   pull(intersection_area)
# 
# # looks good!!
# 
# # # join zcta and state
# # zip_to_state <- st_join(zcta_sf, state_sf, join = st_intersects)
# 
# # join zcta data with patient data
# cco_sf <- right_join(zip_to_state, cco, by = c("ZCTA5CE10" = "zipcode"))
# 
# # get data to plot
# # mean case day temp in each state
# map_temp_df <- cco_sf %>%
#   filter(cases == 1) %>%
#   #filter(STATEFP %in% c("05", "29")) %>%
#   group_by(STATEFP) %>%
#   summarise(mean_temp_lag0 = mean(temp_lag0))
# 
# # temp map
# map_temp_df %>%
#   #filter(STATEFP %in% c("05", "29")) %>%
#   #filter(str_starts(ZCTA5CE10, "0")) %>%
#   #mutate(hosp_here = ifelse()) %>%
#   ggplot() +
#   #geom_sf(aes(fill = temp_lag0), col = NA) +
#   geom_sf(aes(fill = mean_temp_lag0), col = NA) +
#   #geom_sf(data = state_sf, fill = NA) +
#   coord_sf(crs = st_crs(5070)) +
#   theme_void()
# 
# 
# # number of people with hosp
# cco_sf %>%
#   filter(cases == 1) %>%
#   filter(STATEFP %in% c("05", "29")) %>%
#   #filter(str_starts(ZCTA5CE10, "0")) %>%
#   #mutate(hosp_here = ifelse()) %>%
#   ggplot() +
#   #geom_sf(aes(fill = temp_lag0), col = NA) +
#   geom_sf(fill = "blue", col = NA) +
#   #geom_sf(data = state_sf, fill = NA) +
#   coord_sf(crs = st_crs(5070)) +
#   theme_void()
# 
# # now try to figure out what is going on
# # two different geometries? which will it plot?
# # and probably want to do county level, not state



##############################################################

# dataset before trimming 
load(paste0("data/intermediate/cco_", dataset, ".RData"))

# get 95th percentile of 3-day PM2.5
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95) # 18.4 ug/m3

# get dist of PM2.5 before trimming
pdf("results/figures/pm25_dist_pretrim.pdf", width = 8, height = 3)
cco |>
  filter(cases == 1) |>
  ggplot() +
  geom_histogram(aes(pm25_3day),
                 fill = "#c0aada", 
                 binwidth = 0.5) +
  geom_vline(xintercept = pm25_3day_q95, col = "black") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150)) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = expression("Three-day " * PM[2.5] * " (\u03BCg/" * m^3 * ")"),
       y = "Frequency") +
  theme_minimal()
dev.off()

# get range of PM2.5 values (before trimming) on case dates
cases <- cco %>%
  filter(cases == 1)

max(cases$pm25_3day)


