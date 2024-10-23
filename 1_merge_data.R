###############################################################
# Heat-related outcomes
# Authors: Kevin Josey and Lauren Mock
# Aim: Merge data and prepare for analysis
# Inputs: 
#   Exposure data (temp, PM2.5, medication)
#   Outcome data (heat-related hospitalization)
# Output:
#   Case-crossover data for analysis
###############################################################

require(magrittr)
require(haven)
library(readr)
require(lubridate)
require(dplyr)
library(fst)
library(data.table)
library(survival)
#library(dlnm)
library(mgcv)
library(ggplot2)

# choose temperature metric
dataset <- "tmax_temp"
#dataset <- "tmax_pcnt"
#dataset <- "hi_temp"
#dataset <- "hi_pcnt"

#---------- Collect PM2.5 into a single DF ----------#
# Note: don't need to re-run this section (already saved in raw data folder)

# setwd(paste0(proj_path, "data/raw/daily_pm25_sep_files/pm25_all_CSVs"))
# 
# pm_list <- list.files(pattern = ".csv")
# pm_df <- data.frame()
# 
# for (i in 1:length(pm_list)) {
#   
#   pm_tmp <- read_csv(pm_list[[i]])
#   pm_tmp$date <- ymd(substr(pm_list[[i]], 1, 8))
#   pm_df <- rbind(pm_df, pm_tmp)
#   rm(pm_tmp)
#   
# }
# 
# colnames(pm_df) <- tolower(colnames(pm_df))
#
# # Write csv for all PM2.5 data 
# write_csv(pm_df, "data/raw/daily_pm25_one_file.csv")


#---------- Load all data ----------#

##### Temp data

# Read data
temp_df <- read_sas(paste0("data/raw/dlnm_heatrel_", dataset, ".sas7bdat"))


# Change column names
names(temp_df)[3:ncol(temp_df)] <- paste0("temp_", names(temp_df)[3:ncol(temp_df)])
colnames(temp_df) <- tolower(colnames(temp_df))

# Convert to data.table
setDT(temp_df)


##### PM2.5 data

# Read data
pm_df <- read_csv("data/raw/daily_pm25_one_file.csv")

# Change column names
colnames(pm_df) <- tolower(colnames(pm_df))

# Select cols of interest and change names
pm_df <- select(pm_df, -state)
setnames(pm_df, c("zipcode", "pm25_lag0", "dates"))

# Convert to data.table
setDT(pm_df)

# Get PM2.5 lags
setkey(pm_df, zipcode, dates)
pm_df[, paste0("pm25_lag", 1:28) := shift(.SD, 1:28, type = "lag"),
      by = "zipcode", .SDcols = "pm25_lag0"]


##### Medication data

# Binary indicators for medication use each day
# For each individual, we have medication exposures going back 28 days

# Read data
med_df <- read_sas("data/raw/dlnm_med_heatrel.sas7bdat")

# Convert to data.table
setDT(med_df)


##### heat-related event (outcome) data

# Read data
hosp_df <- read_sas("data/raw/cc_heatrel_bi_28d.sas7bdat")

# Change column names
colnames(hosp_df) <- tolower(colnames(hosp_df))

# Factor dow, month, and year
hosp_df$caseyr <- factor(year(hosp_df$dates))
hosp_df$month <- factor(month(hosp_df$dates))
hosp_df$dow <- factor(weekdays(hosp_df$dates), 
                      levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                 "Friday", "Saturday", "Sunday"))

#-- restrict to first hospitalization per beneficiary

# make sure data is ordered properly (by individual and then by da)
hosp_df <- hosp_df[order(hosp_df$bene_id, hosp_df$dates), ]

# get first keyid for each person
first_hosp_keyid <- hosp_df %>%
  group_by(bene_id) %>%
  summarise(first_keyid = first(keyid)) %>%
  pull(first_keyid)

# now filter to rows we care about
hosp_df <- hosp_df %>%
  filter(keyid %in% first_hosp_keyid)


#---------- Set up case-crossover with heat-related hosp data ----------#

# Subset cases
cases <- hosp_df[hosp_df$cases == 1,] 


### Find control days for case days

# Loop through 4 weeks (at most 4 previous Tuesdays in a month)
for (i in 1:4) {
  
  # "Backwards" controls
  # Keep rows where that row's date is i week(s) before the case date and in the same month
  tmp_back <- hosp_df[(hosp_df$dates == hosp_df$casedate - 7*i) &
                        (month(hosp_df$dates) == month(hosp_df$casedate)),]
  
  # "Forwards" controls
  # Keep rows where that row's date is i week(s) after the case date and in the same month
  tmp_for <- hosp_df[(hosp_df$dates == hosp_df$casedate + 7*i) &
                       (month(hosp_df$dates) == month(hosp_df$casedate)),]
  
  # Row bind backwards and forwards controls
  tmp <- bind_rows(tmp_back, tmp_for)
  
  # Row bind all controls from the for loop
  if (i == 1){
    controls <- tmp
  }else{
    controls <- bind_rows(controls, tmp)
  }
  rm(tmp_back, tmp_for, tmp); gc()
}

# Combine cases and controls
cco <- bind_rows(cases, controls)
rm(hosp_df, cases, controls); gc()


# Convert to data.table
setDT(cco)


#---------- Merge with exposures (temp, PM2.5, and steroid use) ----------#

# Merge with temperature
cco <- temp_df[cco, on = .(keyid, dates)]
rm(temp_df); gc()

# Merge with PM2.5
cco <- pm_df[cco, on = .(zipcode, dates)]
rm(pm_df); gc()

# Merge with medication use
cco <- med_df[cco, on = .(keyid)]
rm(med_df); gc()


#---------- Missing data ----------#

# What % of people are missing PM2.5 on lag0?
nrow(cco[cases == 1 & is.na(pm25_lag0)]) / nrow(cco[cases == 1]) 
# 20% missing -- no PM2.5 data for 2017, 2018, 2019

# What % of rows are missing any data? (exclude columns expected to be missing)
cols_to_check <- setdiff(names(cco), c("inhosp", "deathdate",
                                       "med_ANTICHOL_dt", "daysupp_ANTICHOL", "med_STIM_dt",
                                       "daysupp_STIM", "med_LOOPD_dt", "daysupp_LOOPD", "pdx", "sdx",
                                       "user_Antichol", "user_Stim", "user_LoopD"))
nrow(cco[cases == 1 & !complete.cases(cco[, ..cols_to_check]),]) / nrow(cco[cases == 1])

# Same number--only missing data is PM2.5

# Remove these rows with missing data
cco <- cco[complete.cases(cco[, ..cols_to_check]),]


#---------- New columns ----------#

# Column for age (on day of case)
cco[, age := time_length(difftime(casedate, dob), "years") |> trunc()]

# new column for 3-day PM
cco[, pm25_3day := rowMeans(.SD), .SDcols = paste0("pm25_lag", 0:2)]

# column for number of control days per person
cco[, num_control := .N - 1, by = bene_id]


#---------- Save full dataset ----------#

# Save data
save(cco, file = paste0("data/intermediate/cco_", dataset, ".RData"))


#---------- Trim high PM2.5 values and save new dataset ----------#

# get 95th percentile of 3-day PM2.5
pm25_3day_q95 <- quantile(cco$pm25_3day, 0.95) # 18.4 ug/m3

# remove any rows with 3-day PM2.5 > 95th percentile
cco <- cco[pm25_3day <= pm25_3day_q95,]

# remove anyone in the trimmed dataset without a case date
cco <- cco[bene_id %in% cco[cases == 1, unique(bene_id)]]

# new column for number of control days per person after trimming
cco[, num_control_trim := .N - 1, by = bene_id]

# remove anyone with 0 control dates
cco <- cco[num_control_trim != 0,]

# Save trimmed data
save(cco, file = paste0("data/intermediate/cco_trimmed_", dataset, ".RData"))

