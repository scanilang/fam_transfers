library(tidyverse)
library(stats)
library(stargazer)

#######################################################################################
# Read in clean data
#######################################################################################

psid_clean = read.csv('../data/psid_clean.csv')

psid_income_data = psid_clean %>% 
  filter(Age >= 18 & Age <= 60) %>% 
  filter(log_labor_uiwc > 0) 

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
