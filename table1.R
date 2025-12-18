###############################################################
# Temperature, PM2.5, and heat-related hospitalization
# Author: Lauren Mock
# Aim: Create Table 1
# Inputs: 
#   Case-crossover data
# Output:
#   Table 1
###############################################################

library(data.table)
library(table1)
library(tidyverse)
library(xtable)


#---------- load case-crossover data (main analysis) ----------#

# Load case-crossover data
load(paste0("data/intermediate/cco_trim95_tmax_temp.RData"))
cco <- cco_trim95
rm(cco_trim95)


#---------- Convert temperature to celsius ----------#

cco[, temp_lag0 := (temp_lag0 - 32) * (5/9)]


#---------- set up data ----------#

# variables to examine in Table 1
group_vars <- c("sex", "age_cat", "race", "dualeligible")

spacing <- "\\hspace{10pt}"

# new column for age category
cco$age_cat <- cut(cco$age, 
                        breaks = c(65, 74, 84, Inf),
                        labels = paste0(spacing, c("65-74", "75-84", "85+")),
                        right = FALSE)

# change sex, race, dual values
cco[, sex := factor(sex, levels = c(1,2), 
                         labels = paste0(spacing, c("Male", "Female")))]
cco[, race := factor(race, levels = c(1:6), 
                          labels = paste0(spacing, 
                                          c("White", 
                                            "Black/African-American", 
                                            "Other", 
                                            "Asian/Pacific Islander", 
                                            "Hispanic", 
                                            "American Indian/Alaska Native")))]
cco[, dualeligible := factor(dualeligible, levels = c(0,1), 
                                  labels = paste0(spacing, c("Ineligible", "Eligible")))]

# save cases only
cco_cases <- cco[cases == 1]


#---------- beneficiary characteristics ----------#

# write function to get # of people (%) for a bunch of variables (case rows only)
get_n_pcnt <- function(group_vars){
  
  # initialize df
  all_vars_df <- data.frame()
  
  # loop through variables
  for(i in 1:length(group_vars)){
    
    # select column for iteration i
    col <- cco_cases[[group_vars[i]]]
    
    # plot sub-header
    header <- data.frame(var = group_vars[i],
                         group = group_vars[i],
                         value = NA)
    
    # get n for each cat
    one_var_df <- tibble(group = col) %>%
      count(group) %>%
      mutate(pcnt = n / sum(n) * 100,
             var = group_vars[[i]],
             n = comma(n),
             value = paste0(n, " (", round(pcnt, 1), "\\%)")) %>%
      select(var, group, value)
    
    # paste mini tables for all variables together
    all_vars_df <- rbind(all_vars_df, header, one_var_df)
    
  }
  
  return(all_vars_df)
}

# use function
tab1_n_pcnt <- get_n_pcnt(group_vars)

# recode subheaders (now column names)
tab1_n_pcnt <- tab1_n_pcnt %>%
  mutate(group = recode(group, 
                        "age_cat" = "Age category",
                        "race" = "Race/ethnicity",
                        "dualeligible" = "Medicaid eligibility",
                        "sex" = "Sex")) %>%
  select(-var)


#---------- exposures ----------#

# get mean (SD) temperature and PM exposures for case and control days

tab1_expos <- cco %>%
  as.data.frame() %>%
  group_by(cases) %>%
  summarise(temp = paste0(round(mean(temp_lag0), 1), " (", 
                          round(sd(temp_lag0), 1), ")"),
            pm25 = paste0(round(mean(pm25_3day), 1), " (", 
                          round(sd(pm25_3day), 1), ")")) %>%
  pivot_longer(cols = c(temp, pm25),
               names_to = "group") %>%
  mutate(group = paste0(group, "_", cases)) %>%
  select(-cases) %>%
  arrange(desc(group))


#---------- merge and finalize table ----------#

# bind two tables
tab1 <- rbind(tab1_n_pcnt, tab1_expos)

# replace NAs with blanks
tab1 <- tab1 %>%
  mutate_all(~replace_na(., ""))

# get a row for overall
overall_row <- data.frame(
  group = "Total",
  value = format(nrow(cco_cases), big.mark = ","))

# bind to other rows
tab1 <- rbind(overall_row, tab1)

# now clean up complete table 1
tab1 <- tab1 %>%
  rename(" " = group,
         "N (\\%)" = value)

# write to latex
print(xtable(tab1, 
             type = latex,
             label = "tab:table1"), 
      file = "results/tables/table1.tex",
      sanitize.text.function = identity,
      include.rownames = FALSE)
