library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

#######################################################################################
# Read in clean data
#######################################################################################

psid_model_data = read.csv('../../data/psid_model_data.csv')

#######################################################################################
# Step 1: Probit
#######################################################################################

############### separate asset and nonasset income
black_probit_out = glm(Provided_Support_Indicator ~ log_nonasset_income + log_asset_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"),  method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))

black_probit_in = glm(Received_Support_Indicator ~ log_nonasset_income + log_asset_income  + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"),  method = "brglmFit",
                      data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))

white_probit_in = glm(Received_Support_Indicator ~ log_nonasset_income + log_asset_income  + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"), method = "brglmFit",
                      data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))

white_probit_out = glm(Provided_Support_Indicator ~log_nonasset_income + log_asset_income +midpoint_age +  + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"), method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))

stargazer(white_probit_out, black_probit_out,white_probit_in, black_probit_in )

# marginal effects
models <- list(  white_probit_out = white_probit_out,  black_probit_out = black_probit_out, white_probit_in = white_probit_in, black_probit_in = black_probit_in)

mfx_results <- lapply(models, function(model) summary(margins(model)))
options(scipen = 999)
mfx_results_afe = data.frame(mfx_results$white_probit_out$AME,  mfx_results$black_probit_out$AME,
                             mfx_results$white_probit_in$AME,  mfx_results$black_probit_in$AME)

library(xtable)
print(xtable(mfx_results_afe, digits = c(0,4,4,4,4)))

mfx_results_se = data.frame(mfx_results$white_probit_out$SE,  mfx_results$black_probit_out$SE,
                             mfx_results$white_probit_in$SE,  mfx_results$black_probit_in$SE)
print(xtable(mfx_results_se, digits = c(0,4,4,4,4)))
################## asset indicator

black_probit_out = glm(Provided_Support_Indicator ~ log_total_income + as.factor(has_asset_income_3000) + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"),  method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))


black_probit_out = glm(Provided_Support_Indicator ~ log_total_income + log_asset_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"),  method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))

black_probit_in = glm(Received_Support_Indicator ~ log_total_income + as.factor(has_asset_income) + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"),  method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))



summary(black_probit_in)
mfx <- margins(black_probit_in)
summary(mfx)

white_probit_out = glm(Provided_Support_Indicator ~ log_total_income + as.factor(has_asset_income_3000) +midpoint_age +  avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"), method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))



summary(white_probit_out)
mfx <- margins(white_probit_out)
summary(mfx)

stargazer::stargazer(white_probit_out)

white_probit_in = glm(Received_Support_Indicator ~ log_total_income + as.factor(has_asset_income) + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"), method = "brglmFit",
                      data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))




summary(white_probit_in)
mfx <- margins(white_probit_in)
summary(mfx)

# total income only 
############### separate asset and nonasset income
black_probit_out = glm(Provided_Support_Indicator ~ log_total_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"),  method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))

black_probit_in = glm(Received_Support_Indicator ~ log_total_income   + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"),  method = "brglmFit",
                      data = psid_model_data %>% filter(Race_Head == "Black")%>% filter(!is.na(provided_transfer_past)))

white_probit_in = glm(Received_Support_Indicator ~ log_total_income   + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      family = binomial(link = "probit"), method = "brglmFit",
                      data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))

white_probit_out = glm(Provided_Support_Indicator ~log_total_income  +midpoint_age +  + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       family = binomial(link = "probit"), method = "brglmFit",
                       data = psid_model_data %>% filter(Race_Head == "White")%>% filter(!is.na(provided_transfer_past)))

stargazer(white_probit_out, black_probit_out,white_probit_in, black_probit_in )


#######################################################################################
# Step 2: OLS
#######################################################################################

# separate asset and nonasset income
black_transfer_out = lm(log(sum_transfer_out) ~ log_nonasset_income + log_asset_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       data = psid_model_data %>% filter(Race_Head == "Black") %>%  filter(Provided_Support_Indicator == 1))

black_transfer_in = lm(log(sum_transfer_in) ~ log_nonasset_income + log_asset_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      data = psid_model_data %>% filter(Race_Head == "Black") %>%  filter(Received_Support_Indicator == 1))

white_transfer_out = lm(log(sum_transfer_out) ~ log_nonasset_income + log_asset_income + midpoint_age +midpoint_age2 +  avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                       data = psid_model_data %>% filter(Race_Head == "White") %>%  filter(Provided_Support_Indicator == 1))

white_transfer_in = lm(log(sum_transfer_in) ~ log_nonasset_income + log_asset_income + midpoint_age + midpoint_age2 + avg_FU_size + Marital_Status + received_transfer_past + provided_transfer_past + Birth_Cohort + Year_bins,
                      data = psid_model_data %>% filter(Race_Head == "White") %>%  filter(Received_Support_Indicator == 1))

stargazer(white_transfer_out, black_transfer_out,white_transfer_in, black_transfer_in )


