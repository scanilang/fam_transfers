library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

##############################################################################################################################################################################
# Read in clean data
##############################################################################################################################################################################

psid_clean = read.csv('../data/psid_clean.csv') %>% 
  filter(Marital_Status %in% c("Married", "Single"), Age <= 85, Age >= 18)
##############################################################################################################################################################################
# Step 1: Probit
##############################################################################################################################################################################

black_probit_out = glm(Provided_Support_Indicator ~ log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size + 
                         Marital_Status + family_type + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"), 
                       data = psid_clean %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))
summary(black_probit_out)

black_probit_in = glm(Received_Support_Indicator ~ log_nonasset_income + log_asset_income  + Age + Age2 + Family_Unit_Size + 
                        Marital_Status + family_type + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"),  
                      data = psid_clean %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))
summary(black_probit_in)

white_probit_in = glm(Received_Support_Indicator ~ log_nonasset_income + log_asset_income  + Age + Age2 + Family_Unit_Size + 
                        Marital_Status + family_type + received_transfer_past + provided_transfer_past + Birth_Cohort +  Year_bins,
                      family = binomial(link = "probit"), 
                      data = psid_clean %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))
summary(white_probit_in)

white_probit_out = glm(Provided_Support_Indicator ~log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size +
                         Marital_Status + family_type + received_transfer_past + provided_transfer_past + Birth_Cohort +  Year_bins,
                       family = binomial(link = "probit"), 
                       data = psid_clean %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))
summary(white_probit_out)

stargazer(white_probit_out, black_probit_out,white_probit_in, black_probit_in )
# 
# # marginal effects
# models <- list(  white_probit_out = white_probit_out,  black_probit_out = black_probit_out, white_probit_in = white_probit_in, black_probit_in = black_probit_in)
# 
# mfx_results <- lapply(models, function(model) summary(margins(model)))
# options(scipen = 999)
# mfx_results_afe = data.frame(mfx_results$white_probit_out$AME,  mfx_results$black_probit_out$AME,
#                              mfx_results$white_probit_in$AME,  mfx_results$black_probit_in$AME)
# 
# library(xtable)
# print(xtable(mfx_results_afe, digits = c(0,4,4,4,4)))
# 
# mfx_results_se = data.frame(mfx_results$white_probit_out$SE,  mfx_results$black_probit_out$SE,
#                              mfx_results$white_probit_in$SE,  mfx_results$black_probit_in$SE)
# print(xtable(mfx_results_se, digits = c(0,4,4,4,4)))

##############################################################################################################################################################################
# Step 2: OLS
##############################################################################################################################################################################

# separate asset and nonasset income
black_transfer_out = lm(log_transfer_out ~ log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size + 
                          Marital_Status + family_type + Year_bins,
                       data = psid_clean %>% filter(Race_Head == "Black") %>%  filter(Provided_Support_Indicator == 1))
summary(black_transfer_out)

black_transfer_in = lm(log_transfer_in ~ log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size + 
                         Marital_Status + family_type  + Year_bins,
                      data = psid_clean %>% filter(Race_Head == "Black") %>%  filter(Received_Support_Indicator == 1))
summary(black_transfer_in)

white_transfer_out = lm(log_transfer_out ~ log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size + 
                          Marital_Status + family_type  + Year_bins,
                       data = psid_clean %>% filter(Race_Head == "White") %>%  filter(Provided_Support_Indicator == 1))
summary(white_transfer_out)

white_transfer_in = lm(log_transfer_in ~ log_nonasset_income + log_asset_income + Age + Age2 + Family_Unit_Size + 
                         Marital_Status + family_type + Year_bins,
                      data = psid_clean %>% filter(Race_Head == "White") %>%  filter(Received_Support_Indicator == 1))
summary(white_transfer_in)

stargazer(white_transfer_out, black_transfer_out, white_transfer_in, black_transfer_in,
          type = "text" )

#######################################################################################
# Export Results
#######################################################################################


transfer_probit_results = data.frame(white_probit_out = white_probit_out$coefficients,
                                       black_probit_out = black_probit_out$coefficients,
                                       white_probit_in = white_probit_in$coefficients,
                                       black_probit_in = black_probit_in$coefficients)

transfer_amount_results = data.frame(white_transfer_out = white_transfer_out$coefficients,
                                     black_transfer_out = black_transfer_out$coefficients,
                                     white_transfer_in = white_transfer_in$coefficients,
                                     black_transfer_in = black_transfer_in$coefficients)

write.csv(transfer_probit_results, "../data/transfer_probit_results.csv")
write.csv(transfer_amount_results, "../data/transfer_amount_results.csv")


