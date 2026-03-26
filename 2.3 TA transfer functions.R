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

white_probit_edu = glm(Help_Tuition ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type,
                       family = binomial(link = "probit"), 
                       data = tas_clean %>% filter(Race_Head == "White"))
summary(white_probit_edu)

black_probit_edu = glm(Help_Tuition ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type,
                       family = binomial(link = "probit"), 
                       data = tas_clean %>% filter(Race_Head == "Black"))
summary(black_probit_edu)


##############################################################################################################################################################################
# Step 2: OLS
##############################################################################################################################################################################

white_transfer_edu = glm(Help_Tuition_Amount_Parents ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type + degree_type,
                       family = binomial(link = "probit"), 
                       data = tas_clean %>% filter(Race_Head == "White"))
summary(white_transfer_edu)

black_transfert_edu = glm(Help_Tuition_Amount_Parents ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type + degree_type,
                       family = binomial(link = "probit"), 
                       data = tas_clean %>% filter(Race_Head == "Black"))
summary(black_transfert_edu)



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

