library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

##############################################################################################################################################################################
# Read in clean data
##############################################################################################################################################################################

tas_clean = read.csv('../data/tas_clean.csv') 

##############################################################################################################################################################################
# Step 1: Probit
##############################################################################################################################################################################

student_probit = glm(
  Received_Any_Support ~ 
    log_parental_income +
    log_parental_asset_income +
    Race +
    degree_type +
    Cohort,
  family = binomial(link = "probit"),
  data = tas_clean %>% 
    filter(Race %in% c(1,2), 
           Enrollment_Status >= 3,
           Enrollment_Status < 99),
  weights = Individual_Weight
)

##############################################################################################################################################################################
# Step 2: OLS
##############################################################################################################################################################################



#######################################################################################
# Export Results
#######################################################################################


edu_transfer_probit_results = data.frame(white_probit_out = white_probit_out$coefficients,
                                     black_probit_out = black_probit_out$coefficients,
                                     white_probit_in = white_probit_in$coefficients,
                                     black_probit_in = black_probit_in$coefficients)

edu_transfer_amount_results = data.frame(white_transfer_out = white_transfer_out$coefficients,
                                     black_transfer_out = black_transfer_out$coefficients,
                                     white_transfer_in = white_transfer_in$coefficients,
                                     black_transfer_in = black_transfer_in$coefficients)

write.csv(edu_transfer_probit_results, "../data/edu_transfer_probit_results.csv")
write.csv(edu_transfer_amount_results, "../data/edu_transfer_amount_results.csv")

