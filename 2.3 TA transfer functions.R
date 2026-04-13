library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

# Support given when a student

##############################################################################################################################################################################
# Read in clean data
##############################################################################################################################################################################

tas_school_transfers = read.csv('../data/tas_school_transfers.csv') 

##############################################################################################################################################################################
# Step 1: Probit
##############################################################################################################################################################################

white_probit_student_support = glm(Living_Support_Indicator ~ log_nonasset_income + log_asset_income + Head_College,
                       family = binomial(link = "probit"), 
                       data = tas_school_transfers)
summary(white_probit_student_support)

black_probit_student_support = glm(Living_Support_Indicator ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                                     Marital_Status_Parents + family_type,
                       family = binomial(link = "probit"), 
                       data = tas_school_transfers %>% filter(Race_Head == "Black"))
summary(black_probit_student_support)

student_support_probit = glm(Living_Support_Indicator ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                               Marital_Status_Parents + Race_Head,
                             family = binomial(link = "probit"), 
                             data = tas_school_transfers)
summary(student_support_probit)


##############################################################################################################################################################################
# Step 2: OLS
##############################################################################################################################################################################

white_transfer_edu = lm(log(Help_Tuition_Amount_Parents) ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
                         Marital_Status + family_type + degree_type + as.factor(Year),
                       data = tas_clean %>% filter(Race_Head == "White", Help_Tuition_Amount_Parents > 0, !is.na(degree_type)))
summary(white_transfer_edu)

pooled_transfer_edu = lm(
  log(Help_Tuition_Amount_Parents) ~
    log_nonasset_income +
    log_asset_income +
    Head_College +
    family_type +
    degree_type +
    Race_Head +
    as.factor(Year),#+        # identifies intercept shift
 #   Year_bins_edu,
  data = tas_clean %>%
    filter(Race_Head %in% c("White", "Black"),
           Help_Tuition_Amount_Parents > 0,
           !is.na(degree_type),
           !is.na(family_type),
           !is.na(Head_College),
           !is.na(Race_Head),
           !is.na(Year_bins_edu),
           !is.na(Individual_Weight)),
  weights = Individual_Weight
)

summary(pooled_transfer_edu)

# extract race intercept shift
race_shift = coef(pooled_transfer_edu)["Race_HeadWhite"]


# 
# black_transfer_edu = lm(log(Help_Tuition_Amount_Parents) ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College  + degree_type  + as.factor(Year),
#                        data = tas_clean %>% filter(Race_Head == "Black", Help_Tuition_Amount_Parents > 0, !is.na(degree_type)))
# summary(black_transfer_edu)



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

