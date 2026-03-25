library(tidyverse)
library(stats)
library(stargazer)

#######################################################################################
# Read in clean data
#######################################################################################

psid_model_data = read.csv('../../data/psid_model_data.csv')

psid_income_data = psid_model_data %>% 
  filter(midpoint_age >= 22 & midpoint_age <= 60) %>% 
  group_by(Race_Head, Year_first) %>% 
  filter(log_labor_uiwc_income > 0) 

#######################################################################################
# Overlapping panel (STY 2004)
#######################################################################################

# three-period overlapping panel 
psid_income_panel = psid_income_data %>%
  mutate(Panel_Year = Year_first) %>%
  full_join(psid_income_data %>% mutate(Panel_Year = Year_first - 2)) %>%
  full_join(psid_income_data %>% mutate(Panel_Year = Year_first - 4)) %>%
  group_by(ER30001, ER30002, Panel_Year) %>%
  mutate(n = n(),
         lag_income = if_else(lag(Year_first) == Year_first - 2, lag(sum_labor_uiwc), NA),
         growth = (sum_labor_uiwc)/lag_income,
         tag_growth = if_else(abs(growth) > 20 | abs(growth) < (1/20), 1, 0),
         tag_remove = if_else(sum(tag_growth, na.rm  = T) > 0, 1, 0)) %>% 
  ungroup() %>% 
  select(ER30001, ER30002, Panel_Year, n, everything()) %>%
  filter(n == 3, tag_remove == 0) %>%
  arrange(ER30001, ER30002, Panel_Year, Year_first) 


#Each panel begins in a year and consists of observations over that year and the next two years
#when income data are present in all three years.

#######################################################################################
# Deterministic component
#######################################################################################

# PSID data
white_income_process <- lm(log_labor_uiwc_income~ as.factor(Year_first) + midpoint_age + midpoint_age2 + midpoint_age3 + Head_College + avg_FU_size, 
                           data = psid_income_data %>% filter(Race_Head == "White"))

black_income_process <- lm(log_labor_uiwc_income ~ as.factor(Year_first) + midpoint_age + midpoint_age2 + midpoint_age3 + Head_College + avg_FU_size, 
                           data = psid_income_data %>% filter(Race_Head == "Black"))

stargazer(white_income_process, black_income_process)

# Overlapping three-period panel
white_income_process_panel <- lm(log_labor_uiwc_income~ as.factor(Year_first) + midpoint_age + midpoint_age2 + midpoint_age3 + Head_College + avg_FU_size, 
                                 data = psid_income_panel %>% filter(Race_Head == "White"))

black_income_process_panel <- lm(log_labor_uiwc_income~ as.factor(Year_first) + midpoint_age + midpoint_age2 + midpoint_age3 + Head_College + avg_FU_size, 
                                 data = psid_income_panel %>% filter(Race_Head == "Black"))

stargazer(white_income_process_panel, black_income_process_panel)

#######################################################################################
# Idiosyncratic component panel
#######################################################################################

# Year coefficients
year_coefs<-  data.frame(coef= attr(white_income_process_panel$coefficients, "names"), 
                              White = white_income_process_panel$coefficients, 
                              Black = black_income_process_panel$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year_first = if_else(grepl("Year", coef),str_sub(coef, -4,-1),coef ),
         Year_first = if_else(Year_first == "(Intercept)", "1984", Year_first),
         White = if_else(Year_first != "1984", White + white_income_process_panel$coefficients[1], White),
         Black = if_else(Year_first != "1984", Black + black_income_process_panel$coefficients[1], Black)) %>% 
  select(-coef) %>% 
  pivot_longer(!Year_first, names_to = "Race_Head", values_to = "Year_coef") %>% 
  mutate(Year_coef = round(Year_coef, 4),
         Year_first = as.numeric(Year_first)) %>% 
  as.data.frame(.) %>% 
  filter(Year_first %in% 1984:2019)

psid_income_panel = psid_income_panel %>% 
  left_join(year_coefs) %>% # Year coef includes intercept
  mutate(college_coef = case_when(Race_Head == "White" & Head_College == "No College" ~ white_income_process_panel$coefficients["Head_CollegeNo College"],
                                  Race_Head == "White" & Head_College == "Some College" ~ white_income_process_panel$coefficients["Head_CollegeSome College"], 
                                  Race_Head == "Black" & Head_College == "No College" ~ black_income_process_panel$coefficients["Head_CollegeNo College"],
                                  Race_Head == "Black" & Head_College == "Some College" ~ black_income_process_panel$coefficients["Head_CollegeSome College"],
                                  TRUE ~ 0),
          g = if_else(Race_Head == "White", Year_coef +midpoint_age*white_income_process_panel$coefficients["midpoint_age"] + midpoint_age2*white_income_process_panel$coefficients["midpoint_age2"] + 
                       midpoint_age3*white_income_process_panel$coefficients["midpoint_age3"]  + college_coef + 
                       avg_FU_size*white_income_process_panel$coefficients["avg_FU_size"],
                       Year_coef +midpoint_age*black_income_process_panel$coefficients["midpoint_age"] + midpoint_age2*black_income_process_panel$coefficients["midpoint_age2"] + 
                       midpoint_age3*black_income_process_panel$coefficients["midpoint_age3"]  + college_coef + 
                       avg_FU_size*black_income_process_panel$coefficients["avg_FU_size"]),
         u = log_labor_uiwc_income -g) %>% 
  group_by(ER30001, ER30002, Panel_Year) %>% 
  mutate(lagu = lag(u)) %>% 
  ungroup()
  
# AR (1) process
(white_ar1 <- lm(u ~ 0 + lagu, data = psid_income_panel %>% filter(Race_Head == "White")))
(pooled_variance <- var(residuals(white_ar1), na.rm = TRUE))
(sd_ar1 = sqrt(pooled_variance))

(black_ar1 <- lm(u ~ 0 + lagu, data = psid_income_panel %>% filter(Race_Head == "Black")))
(pooled_variance <- var(residuals(black_ar1), na.rm = TRUE))
(sd_ar1 = sqrt(pooled_variance))

#######################################################################################
# Idiosyncratic component no panel
#######################################################################################

# Year coefficients
year_coefs<-  data.frame(coef= attr(white_income_process$coefficients, "names"), 
                         White = white_income_process$coefficients, 
                         Black = black_income_process$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year_first = if_else(grepl("Year", coef),str_sub(coef, -4,-1),coef ),
         Year_first = if_else(Year_first == "(Intercept)", "1984", Year_first),
         White = if_else(Year_first != "1984", White + white_income_process$coefficients[1], White),
         Black = if_else(Year_first != "1984", Black + black_income_process$coefficients[1], Black)) %>% 
  select(-coef) %>% 
  pivot_longer(!Year_first, names_to = "Race_Head", values_to = "Year_coef") %>% 
  mutate(Year_coef = round(Year_coef, 4),
         Year_first = as.numeric(Year_first)) %>% 
  as.data.frame(.) %>% 
  filter(Year_first %in% 1984:2019)

psid_income_data = psid_income_data %>% 
  left_join(year_coefs) %>% # Year coef includes intercept
  mutate(college_coef = case_when(Race_Head == "White" & Head_College == "No College" ~ white_income_process$coefficients["Head_CollegeNo College"],
                                  Race_Head == "White" & Head_College == "Some College" ~ white_income_process$coefficients["Head_CollegeSome College"], 
                                  Race_Head == "Black" & Head_College == "No College" ~ black_income_process$coefficients["Head_CollegeNo College"],
                                  Race_Head == "Black" & Head_College == "Some College" ~ black_income_process$coefficients["Head_CollegeSome College"],
                                  TRUE ~ 0),
         g = if_else(Race_Head == "White", Year_coef +midpoint_age*white_income_process$coefficients["midpoint_age"] + midpoint_age2*white_income_process$coefficients["midpoint_age2"] + 
                       midpoint_age3*white_income_process$coefficients["midpoint_age3"]  + college_coef + 
                       avg_FU_size*white_income_process$coefficients["avg_FU_size"],
                     Year_coef +midpoint_age*black_income_process$coefficients["midpoint_age"] + midpoint_age2*black_income_process$coefficients["midpoint_age2"] + 
                       midpoint_age3*black_income_process$coefficients["midpoint_age3"]  + college_coef + 
                       avg_FU_size*black_income_process$coefficients["avg_FU_size"]),
         u = log_labor_uiwc_income -g,
         Year_lag = Year_first - 2) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lagu = lag(u),
         lagu = if_else(lag(Year_first) ==Year_lag, lagu, NA )) %>% 
  ungroup()

# AR (1) process
(white_ar1 <- lm(u ~ 0 + lagu, data = psid_income_data %>% filter(Race_Head == "White")))
(pooled_variance <- var(residuals(white_ar1), na.rm = TRUE))
(sd_ar1 = sqrt(pooled_variance))

(black_ar1 <- lm(u ~ 0 + lagu, data = psid_income_data %>% filter(Race_Head == "Black")))
(pooled_variance <- var(residuals(black_ar1), na.rm = TRUE))
(sd_ar1 = sqrt(pooled_variance))


########################################################################
#
# Denardi (2004)
#
########################################################################
library(car)
library(sandwich)

data_prep = data_4yr_aggregate %>% 
  mutate(mid_age = case_when(age_cat == 1 ~ 26.5,
                             age_cat == 2 ~ 30.5,
                             age_cat == 3 ~ 34.5,
                             age_cat == 4 ~ 38.5,
                             age_cat == 5 ~ 42.5,
                             age_cat == 6 ~ 46.5,
                             age_cat == 7 ~ 50.5,
                             age_cat == 8 ~ 54.5,
                             age_cat == 9 ~ 58.5,
                             age_cat == 10 ~ 62.5,
                             age_cat == 11 ~ 66.5,
                             age_cat == 12 ~ 70.5,
                             age_cat == 13 ~ 74.5,
                             age_cat == 14 ~ 78.5,
                             age_cat == 15 ~ 82.5),
          mid_age_dif_40 = mid_age - 40, 
          mid_age_dif_40_2 = mid_age_dif_40^2,
          mid_age_dif_40_3 = mid_age_dif_40^3,
         mid_age_2 = mid_age^2,
         mid_age_3 = mid_age^3,
         mid_age_4 = mid_age^4) %>% 
  filter(# !is.na(Completed_Education_max), 
        mid_age <= 58.5)

# specfication 3: regressors are year dummies (take value 1 if the year is the first of a period of 5 years) and fourth-order
# age polynomial, in differences from 40, plus intercept for nonwhite
# plus dummies for less than 5 years of schooling, more than 5 and less than 8, more than 8 and less than 11, 
# college, and more than college

########################################
#
# White
#
########################################
library(stringr)
 
data_prep_white = data_prep %>% filter(Race == "White") 

spec3_white <- lm(log(sum_labor_income_hh) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3 + avg_FU_size, 
                  data = data_prep_white  %>% filter(sum_labor_income_hh >= 4 *1500 & sum_labor_income_hh <= 1800000))

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>%  
         #   Year = if_else(Year == "(Intercept)", '1968', Year),
       #  Year = as.integer(Year)) %>% 
  select(Year, year_coef)

data_final_white = data_prep_white %>% 
  filter(sum_labor_income_hh >= 4 *1500) %>% 
  #mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(log_inc = log(sum_labor_income_hh) ,
         g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"]  + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_hh) - g) %>% 
  arrange(ER30001, ER30002, age_cat) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number(),
         count = n(),
         count_lag = sum(!is.na(lag_u))) %>%  
  ungroup() %>% 
  mutate(id = paste0(ER30001, " ", ER30002))


##### ATTEMPT 2
ar1_model <- lm(u ~ 0 + lag_u, data = data_final_white)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))

spec3_white <- lm(log(sum_labor_income_comp_ui) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3 + avg_FU_size, 
                  data = data_prep_white %>% filter(sum_labor_income_comp_ui >= 4 *1500 &sum_labor_income_comp_ui <= 1800000) )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>%  
  #   Year = if_else(Year == "(Intercept)", '1968', Year),
  #  Year = as.integer(Year)) %>% 
  select(Year, year_coef)

data_final_white = data_prep_white %>%  filter(sum_labor_income_comp_ui >= 4 * 1500 &sum_labor_income_comp_ui <= 1500000) %>% 
  #mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(log_inc = log(sum_labor_income_comp_ui) ,
         g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"]  + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_comp_ui) - g,
         log = g + u) %>% 
  arrange(ER30001, ER30002, age_cat) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number(),
         count = n(),
         count_lag = sum(!is.na(lag_u))) %>%  
  ungroup() %>% 
  mutate(id = paste0(ER30001, " ", ER30002))

ar1_model <- lm(u ~ 0 + lag_u, data = data_final_white)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))

spec3_white <- lm(log(sum_labor_income_head) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3 + avg_FU_size, 
                  data = data_prep_white %>% filter(sum_labor_income_head > 4*1500) )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>%  
  #   Year = if_else(Year == "(Intercept)", '1968', Year),
  #  Year = as.integer(Year)) %>% 
  select(Year, year_coef)

data_final_white = data_prep_white %>% 
 filter(sum_labor_income_head > 4*1500) %>% 
  #mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(log_inc = log(sum_labor_income_head) ,
         g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"]  + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_head) - g) %>% 
  arrange(ER30001, ER30002, age_cat) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number(),
         count = n(),
         count_lag = sum(!is.na(lag_u))) %>%  
  ungroup() %>% 
  mutate(id = paste0(ER30001, " ", ER30002))

ar1_model <- lm(u ~ 0 + lag_u, data = data_final_white)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))

## lag 1 autocorrelation
autocor_white <- data_final_white %>% 
  filter(count_lag > 0) %>% 
#  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u),
         u_sub_mean = u - mean_u,
         mean_sub = mean(u_sub_mean),
         lagu_sub_mean = lag_u - mean_u,
         top_element = u_sub_mean * lagu_sub_mean,
         bottom_element = u_sub_mean^2) %>% 
  ungroup() %>% 
  mutate(var_ucenter = var(u_sub_mean),
         mean_ucenter = mean(u_sub_mean)) %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(covariance =  mean(lagu_sub_mean* u_sub_mean, na.rm = TRUE),
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator,
         r2 = covariance / var_u ) 

# moments from data
autocor_white$mean_ucenter[1] # sample mean u should be close to 0
var(autocor_white$u_sub_mean) # sample variance
(var_ucenter  = autocor_white$var_ucenter[1]) # sample variance
(rho_white = autocor_white$r1[1])
(rho_white = autocor_white$r2[1])
(sigma2 = var_ucenter[1]  * (1-(rho_white)^2))


##### USING INDIVIDUAL YEARS

spec3_white <- lm(log(sum_labor_income_hh) ~ as.factor(Year_first) + mid_age + mid_age_2 + mid_age_3 + avg_FU_size, 
                  data = data_prep_white )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -4,-1),coef ),
    Year = if_else(Year == "(Intercept)", '1968', Year),
    Year = as.integer(Year)) %>% 
  select(Year, year_coef)

data_final_white = data_prep_white %>% 
  #mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first' = 'Year')) %>% 
  mutate(log_inc = log(sum_labor_income_hh) ,
         g = if_else(Year_first != "1984", spec3_white$coefficients[1] + year_coef  +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"]  + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_hh) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number(),
         count = n(),
         count_lag = sum(!is.na(lag_u))) %>%  
  ungroup()

## lag 1 autocorrelation
autocor_white <- data_final_white %>% 
  filter(count_lag > 0) %>% 
  #  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u),
         u_sub_mean = u - mean_u,
         mean_sub = mean(u_sub_mean),
         lagu_sub_mean = lag_u - mean_u,
         top_element = u_sub_mean * lagu_sub_mean,
         bottom_element = u_sub_mean^2) %>% 
  ungroup() %>% 
  mutate(var_ucenter = var(u_sub_mean),
         mean_ucenter = mean(u_sub_mean)) %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(covariance =  mean(lagu_sub_mean* u_sub_mean, na.rm = TRUE),
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator,
         r2 = covariance / var_u ) 

# moments from data
autocor_white$mean_ucenter[1] # sample mean u should be close to 0
var(autocor_white$u_sub_mean) # sample variance
(var_ucenter  = autocor_white$var_ucenter[1]) # sample variance
(rho_white = autocor_white$r1[1])
(rho_white = autocor_white$r2[1])
(sigma2 = var_ucenter[1]  * (1-(rho_white)^2))

  
########################################
#
# Black
#
########################################
data_prep_black = data_prep %>% filter(Race == "Black")

spec3_black <- lm(log(sum_labor_income_hh) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3 +avg_FU_size, 
                  data = data_prep_black %>%  filter(sum_labor_income_hh >= 4 *1500 & sum_labor_income_hh <= 1800000))

year_coef_data <-  data.frame(coef= attr(spec3_black$coefficients, "names"), year_coef = spec3_black$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_black = data_prep_black %>%  filter(sum_labor_income_hh >= 4 *1500) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != '1984-1990', spec3_black$coefficients[1] + year_coef  +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"],
                     spec3_black$coefficients[1] +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_hh) - g,
         log = g + u) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number()) %>%  
  ungroup()

ar1_model <- lm(u ~ 0 + lag_u, data = data_final_black)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))


spec3_black <- lm(log(sum_labor_income_head) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3   +avg_FU_size, 
                  data = data_prep_black %>% filter(sum_labor_income_head > 4*1000))

year_coef_data <-  data.frame(coef= attr(spec3_black$coefficients, "names"), year_coef = spec3_black$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_black = data_prep_black %>%
 filter(sum_labor_income_head > 4*1500) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != '1984-1990', spec3_black$coefficients[1] + year_coef  +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"],
                     spec3_black$coefficients[1] +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_head) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number()) %>%  
  ungroup()

ar1_model <- lm(u ~ 0 + lag_u, data = data_final_black)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))


spec3_black <- lm(log(sum_labor_income_comp_ui) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3 +avg_FU_size, 
                  data = data_prep_black %>% filter(sum_labor_income_comp_ui >= 4*1500 & sum_labor_income_comp_ui <= 1800000))

year_coef_data <-  data.frame(coef= attr(spec3_black$coefficients, "names"), year_coef = spec3_black$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_black = data_prep_black %>%
  filter(sum_labor_income_comp_ui >= 4*1500) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != '1984-1990', spec3_black$coefficients[1] + year_coef  +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"],
                     spec3_black$coefficients[1] +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"] + 
                       avg_FU_size*spec3_black$coefficients["avg_FU_size"]),
         u = log(sum_labor_income_comp_ui) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number()) %>%  
  ungroup()

ar1_model <- lm(u ~ 0 + lag_u, data = data_final_black)
summary(ar1_model)
(pooled_variance <- var(residuals(ar1_model), na.rm = TRUE))

#using year bins
spec3_black <- lm(log(sum_fam_income) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3  + avg_, 
                  data = data_prep_black )

year_coef_data <-  data.frame(coef= attr(spec3_black$coefficients, "names"), year_coef = spec3_black$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_black = data_prep_black %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != "1984-1990", spec3_black$coefficients[1] + year_coef  +
                       mid_age*spec3_black$coefficients["mid_age"] + mid_age_2*spec3_black$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_black$coefficients["mid_age_3"]  + #Completed_Education_max*spec3_white$coefficients["Completed_Education_max"]+ 
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + #Completed_Education_max*spec3_white$coefficients["Completed_Education_max"]+
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"]),
         u = log(sum_fam_income) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_age = lag(age_cat),
         lag_u = if_else(age_cat == lag_age + 1, lag(u), NA),
         number = row_number(),
         count = n()) %>%  
  ungroup()

## lag 1 autocorrelation
autocor_black <- data_final_black %>% 
  #filter(count > 1) %>% 
  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u)) %>% 
  ungroup() %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(u_sub_mean = u - mean_u,
         lagu_sub_mean = lag_u - mean_u,
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator) 

# moments from data
mean(data_final_black$u) # sample mean u should be close to 0
var(autocor_black$u) # sample variance
(rho_black = autocor_black$r1[1])
sigma2_black = var(autocor_white$u)  * (1-(rho_black)^2)
(sigma_black = sqrt(sigma2))

write.csv(data_final_black, "data clean/black_income_data.csv")


########################################################################
#
# Individual Years
#
########################################################################
data_prep_individ = data_4yr %>% 
  mutate(log_total_fam = log(Total_Family_Income),
    mid_age = case_when(age_cat == 1 ~ 26.5,
                        age_cat == 2 ~ 30.5,
                        age_cat == 3 ~ 34.5,
                        age_cat == 4 ~ 38.5,
                        age_cat == 5 ~ 42.5,
                        age_cat == 6 ~ 46.5,
                        age_cat == 7 ~ 50.5,
                        age_cat == 8 ~ 54.5,
                        age_cat == 9 ~ 58.5,
                        age_cat == 10 ~ 62.5,
                        age_cat == 11 ~ 66.5,
                        age_cat == 12 ~ 70.5,
                        age_cat == 13 ~ 74.5,
                        age_cat == 14 ~ 78.5,
                        age_cat == 15 ~ 82.5),
    mid_age_dif_40 = mid_age - 40, 
    mid_age_dif_40_2 = mid_age_dif_40^2,
    mid_age_dif_40_3 = mid_age_dif_40^3,
    mid_age_2 = mid_age^2,
    mid_age_3 = mid_age^3,
    mid_age_4 = mid_age^4,
    Age2 = Age^2,
    Age3 = Age^3,
    Age4 = Age^4,
    id = paste0(ER30001, " ", ER30002)) %>% 
  filter(
         !is.na(age_cat),mid_age <= 58.5, Labor_Income_Head > 1000)

### WHITE
spec3_white <- lm(log(Labor_Income_Head) ~ as.factor(Year_first_bin) + mid_age + mid_age_2 + mid_age_3  + mid_age_4 + avg_FU_size, 
                  data = data_prep_individ %>%filter(Race == "White") )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_white = data_prep_individ %>%filter(Race == "White")  %>% 
  mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"]  + mid_age_4*spec3_white$coefficients["mid_age_4"] + 
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"],
                     spec3_white$coefficients[1] +
                       mid_age*spec3_white$coefficients["mid_age"] + mid_age_2*spec3_white$coefficients["mid_age_2"] + 
                       mid_age_3*spec3_white$coefficients["mid_age_3"] + mid_age_4*spec3_white$coefficients["mid_age_4"] +
                       avg_FU_size*spec3_white$coefficients["avg_FU_size"]),
         u = log(Labor_Income_Head) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_year = lag(Year),
         lag_u = if_else(Year == lag_year + 1, lag(u), NA),
         number = row_number(),
         count = n()) %>%  
  ungroup()

## lag 1 autocorrelation
autocor_white <- data_final_white %>% 
  filter(count > 4) %>% 
  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u),
         u_sub_mean = u - mean_u,
         lagu_sub_mean = lag_u - mean_u,
         top_element = u_sub_mean * lagu_sub_mean,
         bottom_element = u_sub_mean^2) %>% 
  ungroup() %>% 
  mutate(var_u = var(u_sub_mean)) %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(covariance =  mean(lagu_sub_mean* u_sub_mean, na.rm = TRUE),
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator,
         r2 = covariance / var_u ) 

# moments from data
mean(data_final_white$u) # sample mean u should be close to 0
var(autocor_white$u_sub_mean) # sample variance
(rho_white = autocor_white$r1[1])
(rho_white = autocor_white$r2[1])
(sigma2 = var(autocor_white$u_sub_mean)  * (1-(rho_white)^2))
(sigma2 = var(autocor_white$u)  * (1-(rho_white)^2))

### USING INDIVIDUAL AGE
spec3_white <- lm(log(Total_Family_Income) ~ as.factor(Year_first_bin) + Age + Age2 + Age3  + Age4 + Family_Unit_Size, 
                  data = data_prep_individ %>%filter(Race == "White") )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_white = data_prep_individ %>%filter(Race == "White")  %>% 
  mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       Age*spec3_white$coefficients["Age"] + Age2*spec3_white$coefficients["Age2"] + 
                       Age3*spec3_white$coefficients["Age3"]  + Age4*spec3_white$coefficients["Age4"] + 
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"],
                     spec3_white$coefficients[1] +
                       Age*spec3_white$coefficients["Age"] + Age2*spec3_white$coefficients["Age2"] + 
                       Age3*spec3_white$coefficients["Age3"] + Age4*spec3_white$coefficients["Age4"] +
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"]),
         u = log(sum_fam_income) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_year = lag(Year),
         lag_u = if_else(Year == lag_year + 1, lag(u), NA),
         number = row_number(),
         count = n()) %>%  
  ungroup()

## lag 1 autocorrelation
autocor_white <- data_final_white %>% 
  filter(count > 4) %>% 
  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u),
         u_sub_mean = u - mean_u,
         lagu_sub_mean = lag_u - mean_u,
         top_element = u_sub_mean * lagu_sub_mean,
         bottom_element = u_sub_mean^2) %>% 
  ungroup() %>% 
  mutate(var_u = var(u_sub_mean)) %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(covariance =  mean(lagu_sub_mean* u_sub_mean, na.rm = TRUE),
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator,
         r2 = covariance / var_u ) 

# moments from data
mean(data_final_white$u) # sample mean u should be close to 0
var(autocor_white$u_sub_mean) # sample variance
(rho_white = autocor_white$r1[1])
(rho_white = autocor_white$r2[1])
(sigma2 = var(autocor_white$u_sub_mean)  * (1-(rho_white)^2))
(sigma = sqrt(sigma2))
(sigma2 = var(autocor_white$u)  * (1-(rho_white)^2))

## Black

spec3_black <- lm(log(Total_Family_Income) ~ as.factor(Year_first_bin) + Age + Age2 + Age3  + Age4 + Family_Unit_Size, 
                  data = data_prep_individ %>%filter(Race == "Black") )

year_coef_data <-  data.frame(coef= attr(spec3_white$coefficients, "names"), year_coef = spec3_white$coefficients) %>%  
  tibble::remove_rownames() %>% 
  mutate(Year = if_else(grepl("Year", coef),str_sub(coef, -9,-1),coef )) %>% 
  select(Year, year_coef)

data_final_black = data_prep_individ %>%filter(Race == "Black")  %>% 
  mutate(Married = if_else(Marital_Status == "Single", 1, 0)) %>% 
  left_join(year_coef_data, by = c('Year_first_bin' = 'Year')) %>% 
  mutate(g = if_else(Year_first_bin != "1984-1990", spec3_white$coefficients[1] + year_coef  +
                       Age*spec3_white$coefficients["Age"] + Age2*spec3_white$coefficients["Age2"] + 
                       Age3*spec3_white$coefficients["Age3"]  + Age4*spec3_white$coefficients["Age4"] + 
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"],
                     spec3_white$coefficients[1] +
                       Age*spec3_white$coefficients["Age"] + Age2*spec3_white$coefficients["Age2"] + 
                       Age3*spec3_white$coefficients["Age3"] + Age4*spec3_white$coefficients["Age4"] +
                       Family_Unit_Size*spec3_white$coefficients["Family_Unit_Size"]),
         u = log(sum_fam_income) - g) %>% 
  group_by(ER30001, ER30002) %>% 
  mutate(lag_year = lag(Year),
         lag_u = if_else(Year == lag_year + 1, lag(u), NA),
         number = row_number(),
         count = n()) %>%  
  ungroup()

## lag 1 autocorrelation
autocor_black <- data_final_black %>% 
  filter(count > 4) %>% 
  filter(!is.na(lag_u) | number == 1 ) %>% ### NEED TO INCLUDE MEAN OF THE FIRST PERIOD
  group_by(ER30001, ER30002) %>% 
  mutate(mean_u = mean(u),
         u_sub_mean = u - mean_u,
         lagu_sub_mean = lag_u - mean_u,
         top_element = u_sub_mean * lagu_sub_mean,
         bottom_element = u_sub_mean^2) %>% 
  ungroup() %>% 
  mutate(var_u = var(u_sub_mean)) %>% 
  filter(!is.na(lag_u)) %>% 
  mutate(covariance =  mean(lagu_sub_mean* u_sub_mean, na.rm = TRUE),
         numerator = sum(u_sub_mean * lagu_sub_mean),
         denominator = sum(u_sub_mean^2),
         r1 = numerator / denominator,
         r2 = covariance / var_u ) 

# moments from data
mean(data_final_black$u) # sample mean u should be close to 0
var(autocor_white$u_sub_mean) # sample variance
var(autocor_white$u) # sample variance
(rho_black = autocor_black$r1[1])
(rho_black = autocor_black$r2[1])
(sigma2 = var(autocor_black$u_sub_mean)  * (1-(rho_black)^2))
(sigma = sqrt(sigma2))

(sigma2 = var(autocor_black$u)  * (1-(rho_black)^2))




##### 

y = mean(data_prep$sum_labor_income_head, na.rm =T)
t2 = 2908
0.258*(y - (y^(-0.768) + t2)^(1/0.768))
57983.2

14495.8


### 
tax_y <- function(y,t){
  val = 0.258*(y - (y^(-0.768) + t)^(1/0.768))
  return(val)
}
  
# tax 2010 bracket
tax_amount = (y - 8400 - 45550) * 0.25 + 6235 # head of household
tax_amount = (y - 34000) * 0.25 + 4681.25 # unmarried
