library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

#################################################################################################################################
# Read in clean data
#################################################################################################################################

tas_clean = read.csv('../data/tas_clean.csv') 
psid_edu = read.csv('../data/psid_edu.csv') 

#################################################################################################################################
# Step 1: Probit
#################################################################################################################################

# using TAS
# white_probit_edu = glm(Help_Tuition ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
#                          Marital_Status + family_type + Year_bins_edu,
#                        family = binomial(link = "probit"), 
#                        data = tas_clean %>% filter(Race_Head == "White", Enrollment_Status %in% c(9, 10),  # currently enrolled
#                                                    Part_or_Full_Time_Student == 1))
# summary(white_probit_edu)
# 
# black_probit_edu = glm(Help_Tuition ~ log_nonasset_income + log_asset_income  + Head_College  + Year_bins_edu,
#                        family = binomial(link = "probit"), 
#                        data = tas_clean %>% filter(Race_Head == "Black", Enrollment_Status %in% c(9, 10),  # currently enrolled
#                                                    Part_or_Full_Time_Student == 1))
# summary(black_probit_edu)


# using family roster and transfers
# problem : i dont know part or full time
white_probit_edu = glm(Help_School_Indicator ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type  ,
                       family = binomial(link = "probit"), 
                       data = psid_edu %>% filter(Race_Head == "White"))
summary(white_probit_edu)

black_probit_edu = glm(Help_School_Indicator ~ log_nonasset_income + log_asset_income  + Head_College  ,
                       family = binomial(link = "probit"), 
                       data = psid_edu %>% filter(Race_Head == "Black"))
summary(black_probit_edu)



#################################################################################################################################
# Step 2: OLS
#################################################################################################################################

# white_transfer_edu <- lm(log_educ_exp ~ log_nonasset_income + log_asset_income +
#                            Family_Unit_Size + Head_College + Marital_Status + family_type +
#                            degree_type +enroll_era, 
#                          data    = psid_edu %>% filter(Race_Head == "White", Help_School_Indicator ==1 ))
# summary(white_transfer_edu)
# 
# black_transfer_edu <- lm(log_educ_exp ~ log_nonasset_income + log_asset_income +
#                            Family_Unit_Size + Head_College + Marital_Status + family_type +
#                            degree_type + enroll_era,
#                          data    = psid_edu %>% filter(Race_Head == "Black", Help_School_Indicator ==1 ))
# summary(black_transfer_edu)

pooled_transfer_edu <- lm(
  log_educ_exp ~ log_nonasset_income + log_asset_income +
    Head_College + degree_type + enroll_era +
    Race_Head,
  data = psid_edu %>% filter(Help_School_Indicator == 1,
                             degree_type %in% c("2yr", "4yr"),
                             Race_Head %in% c("White", "Black")))
summary(pooled_transfer_edu)

#################################################################################################################################
# Export Results
#################################################################################################################################


edu_transfer_probit_results = data.frame(white_probit_edu = white_probit_edu$coefficients,
                                     black_probit_edu = black_probit_edu$coefficients)

edu_transfer_amount_results = data.frame(pooled_transfer_edu = pooled_transfer_edu$coefficients)

write.csv(edu_transfer_probit_results, "../data/edu_transfer_probit_results.csv")
write.csv(edu_transfer_amount_results, "../data/edu_transfer_amount_results.csv")

