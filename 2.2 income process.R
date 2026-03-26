library(tidyverse)
library(stats)
library(stargazer)

#######################################################################################
# Read in clean data
#######################################################################################

psid_clean = read.csv('../data/psid_clean.csv')

psid_income_data = psid_clean %>% 
  filter(Age >= 18 & Age <= 60) %>% 
  filter(log_labor_uiwc > 0)  %>% 
 # filter( !(Employment_Status_Ind %in% 4:7))  %>%   # student, not permanently disabled, retired, house
  mutate(Age3 = Age^3,
         Marital_Status_Ind = if_else(Marital_Status == "Married", 1, 0)) 

#######################################################################################
# Deterministic component
#######################################################################################

# PSID data
white_income_process <- lm(log_labor_uiwc~  Age + Age2 + Age3 + Head_College + Marital_Status  + as.factor(Year), 
                           data = psid_income_data %>% filter(Race_Head == "White", Labor_UIWC_Income > 500))

black_income_process <- lm(log_labor_uiwc ~ Age + Age2 + Age3 + Head_College + Marital_Status  + as.factor(Year), 
                           data = psid_income_data %>% filter(Race_Head == "Black", Labor_UIWC_Income > 500))

summary(white_income_process)
summary(black_income_process)

#######################################################################################
# Idiosyncratic component
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

psid_income_panel = psid_income_data %>% 
  left_join(year_coefs) %>% # Year coef includes intercept
  mutate(college_coef = case_when(Race_Head == "White" & Head_College == "No College" ~ white_income_process$coefficients["Head_CollegeNo College"],
                                  Race_Head == "White" & Head_College == "Some College" ~ white_income_process$coefficients["Head_CollegeSome College"], 
                                  Race_Head == "Black" & Head_College == "No College" ~ black_income_process$coefficients["Head_CollegeNo College"],
                                  Race_Head == "Black" & Head_College == "Some College" ~ black_income_process$coefficients["Head_CollegeSome College"],
                                  TRUE ~ 0),
          g = if_else(Race_Head == "White", Year_coef +Age*white_income_process$coefficients["Age"] + Age2*white_income_process$coefficients["Age2"] + 
                        Age3*white_income_process$coefficients["Age3"]  + college_coef + 
                        Marital_Status_Ind*white_income_process$coefficients["Marital_StatusMarried"],
                       Year_coef +Age*black_income_process$coefficients["Age"] + Age2*black_income_process$coefficients["Age2"] + 
                        Age3*black_income_process$coefficients["Age3"]  + college_coef + 
                        Marital_Status_Ind*black_income_process$coefficients["Marital_StatusMarried"]),
         u = log_labor_uiwc -g) %>% 
  group_by(ER30001, ER30002) %>% 
  arrange(ER30001, ER30002, Survey_Year) %>%
  mutate(lagu = lag(u),
         lag_years = Survey_Year - lag(Survey_Year)) %>% 
  ungroup()
  
# AR (1) process
# 1-year lag (pre-1997 annual waves)
white_ar1_1yr <- lm(u ~ 0 + lagu,
                    data = psid_income_panel %>%
                      filter(Race_Head == "White",
                             lag_years == 1))

# 2-year lag (post-1997 biannual waves)
white_ar1_2yr <- lm(u ~ 0 + lagu,
                    data = psid_income_panel %>%
                      filter(Race_Head == "White",
                             lag_years == 2))

# 1-year lag (pre-1997 annual waves)
black_ar1_1yr <- lm(u ~ 0 + lagu,
                    data = psid_income_panel %>%
                      filter(Race_Head == "Black",
                             lag_years == 1))

# 2-year lag (post-1997 biannual waves)
black_ar1_2yr <- lm(u ~ 0 + lagu,
                    data = psid_income_panel %>%
                      filter(Race_Head == "Black",
                             lag_years == 2))

# annual persistence
rho_white <- sqrt(0.651)   # = 0.807
rho_black <- sqrt(0.596)   # = 0.772


# annual shock variance
var_2yr_white <- var(residuals(white_ar1_2yr), na.rm = TRUE)
var_2yr_black <- var(residuals(black_ar1_2yr), na.rm = TRUE)

sigma_white <- sqrt(var_2yr_white / (1 + rho_white^2))
sigma_black <- sqrt(var_2yr_black / (1 + rho_black^2))

# summary table
data.frame(
  Race  = c("White", "Black"),
  rho   = round(c(rho_white, rho_black), 3),
  sigma = round(c(sigma_white, sigma_black), 3)
)


# (white_ar1 <- lm(u ~ 0 + lagu, data = psid_income_panel %>% filter(Race_Head == "White")))
# (pooled_variance <- var(residuals(white_ar1), na.rm = TRUE))
# (sd_ar1 = sqrt(pooled_variance))
# 
# (black_ar1 <- lm(u ~ 0 + lagu, data = psid_income_panel %>% filter(Race_Head == "Black")))
# (pooled_variance <- var(residuals(black_ar1), na.rm = TRUE))
# (sd_ar1 = sqrt(pooled_variance))


#######################################################################################
# Export Results
#######################################################################################
income_results = data.frame(white_income_process = white_income_process$coefficients,
                            black_income_process = black_income_process$coefficients)

write.csv(income_results, "../data/income_results.csv")


#######################################################################################
# Consumption Floor
#######################################################################################

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
