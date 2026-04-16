library(tidyverse)
library(margins)
library(brglm2)
library(stargazer)

#################################################################################################################################
# Read in clean data
#################################################################################################################################

tas_clean = read.csv('../data/tas_clean.csv') 
psid_edu = read.csv('../data/psid_edu.csv') %>% 
  mutate(Head_College = factor(Head_College, 
                               levels = c("No College", "Some College", "College Degree")),
         Marital_Status = factor(Marital_Status,
                                 levels = c("Single", "Married")),
         family_type = factor(family_type,
                              levels = c("low", "mid", "high")))


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
# white_probit_edu = glm(Help_School_Indicator ~ log_nonasset_income + log_asset_income + Family_Unit_Size + Head_College +
#                          Marital_Status + family_type  ,
#                        family = binomial(link = "probit"), 
#                        data = psid_edu %>% filter(Race_Head == "White", degree_type_final %in% c("2yr", "4yr")))
# summary(white_probit_edu)
# 
# black_probit_edu = glm(Help_School_Indicator ~ log_nonasset_income + log_asset_income  +  Family_Unit_Size + Head_College +
#                          Marital_Status   ,
#                        family = binomial(link = "probit"), 
#                        data = psid_edu %>% filter(Race_Head == "Black", degree_type_final %in% c("2yr", "4yr")))
# summary(black_probit_edu)

pooled_probit_edu = glm(Help_School_Indicator ~ log_nonasset_income + log_asset_income + degree_type_final+
                          Family_Unit_Size + Head_College + Marital_Status + family_type +
                          Race_Head + enroll_era,
                        family = binomial(link = "probit"), 
                        data = psid_edu %>% filter(Race_Head %in% c("White", "Black"), Marital_Status %in% c("Married", "Single"),
                                                   degree_type_final %in% c("2yr", "4yr")))
summary(pooled_probit_edu)


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
    Head_College  + degree_type_final + enroll_era + 
    Race_Head,
  data = psid_edu %>% filter(Help_School_Indicator == 1,
                             degree_type_final %in% c("2yr", "4yr"), Marital_Status %in% c("Married", "Single"),
                             Race_Head %in% c("White", "Black")))
summary(pooled_transfer_edu)

#################################################################################################################################
# Export Results
#################################################################################################################################

# 
# edu_transfer_probit_results = data.frame(white_probit_edu = white_probit_edu$coefficients,
#                                      black_probit_edu = black_probit_edu$coefficients)

edu_transfer_probit_results = data.frame(pooled_probit_edu = pooled_probit_edu$coefficients)
edu_transfer_amount_results = data.frame(pooled_transfer_edu = pooled_transfer_edu$coefficients)

write.csv(edu_transfer_probit_results, "../data/edu_transfer_probit_results.csv")
write.csv(edu_transfer_amount_results, "../data/edu_transfer_amount_results.csv")

stargazer(pooled_probit_edu, pooled_transfer_edu,
          type = "latex",
          title = "Education Transfer Regressions",
          column.labels = c("Probit", "OLS"),
          dep.var.labels = c("P(Transfer)", "Log Amount"),
          covariate.labels = c("Log Non-Asset Income",
                               "Log Asset Income",
                               "Degree Type: 4yr",
                               "Family Unit Size",
                               "Parental Educ: Some Colleg",
                               "Parental Educ: College Degree",
                               "Marital Status: Married",
                               "Family Type: Mid",
                               "Family Type: High",
                               "Race: White",
                               "Enrollment Era: 2009-2012",
                               "Enrollment Era: pre-2002"),
          keep.stat = c("n", "rsq", "adj.rsq", "ll"),
          notes = "Sample: White and Black families, 2yr and 4yr degree types. Amount deflated to 2010\\$.",
          notes.append = FALSE,
          label = "tab:edu_transfers",
          out = "edu_transfer_table.tex")

#################################################################################################################################
# Graphs
#################################################################################################################################

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(Help_School_Indicator)) %>%
  group_by(Race_Head, degree_type_final) %>%
  summarize(share = mean(Help_School_Indicator, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = degree_type_final, y = share, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Share Receiving Education Transfers (RT13)",
       x = "Degree Type", y = "Share Receiving", fill = "Race") +
  theme_minimal()

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         Help_School_Indicator == 1,
         Help_School_Amount_Adj > 0) %>%
  group_by(Race_Head, degree_type_final) %>%
  summarize(mean_amt = mean(Help_School_Amount_Adj, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = degree_type_final, y = mean_amt, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::dollar) +
  labs(title = "Mean Education Transfer Amount, Conditional on Receiving (2010 $)",
       x = "Degree Type", y = "Mean Amount", fill = "Race") +
  theme_minimal()

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(log_nonasset_income)) %>%
  mutate(income_quintile = ntile(log_nonasset_income, 5)) %>%
  group_by(Race_Head, income_quintile) %>%
  summarize(share = mean(Help_School_Indicator, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = income_quintile, y = share, color = Race_Head, group = Race_Head)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Transfer Receipt by Parental Income Quintile",
       x = "Income Quintile (1 = Lowest)", y = "Share Receiving", color = "Race") +
  theme_minimal()

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(Head_College)) %>% 
  mutate(Head_College = factor(Head_College, c("No College", "Some College", "College Degree"))) %>%
  group_by(Race_Head, Head_College) %>%
  summarize(share = mean(Help_School_Indicator, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Head_College, y = share, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Transfer Receipt by Parental Education",
       x = "Head College Status", y = "Share Receiving", fill = "Race") 

ggsave("edu_share_transfer_parent_edu.png")


psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(log_nonasset_income),
         Help_School_Indicator == 1,
         Help_School_Amount_Adj > 0) %>%
  mutate(Head_College = factor(Head_College, c("No College", "Some College", "College Degree"))) %>%
  group_by(Race_Head, Head_College) %>%
  summarize(mean_amt = mean(Help_School_Amount_Adj, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(Head_College), y = mean_amt, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::dollar) +
  labs(title = "Mean Education Transfer Amount by Parental Education",
       subtitle = "Conditional on Receiving a Transfer (2010 $)",
       x = "Head College Status", y = "Mean Amount", fill = "Race") 

ggsave("edu_mean_transfer_parent_edu.png")

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(log_nonasset_income)) %>%
  mutate(income_quintile = ntile(log_nonasset_income, 4)) %>%
  group_by(Race_Head, income_quintile) %>%
  summarize(share = mean(Help_School_Indicator, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(income_quintile), y = share, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Transfer Receipt by Parental Income Quantile",
       x = "Income Quantile (1 = Lowest)", y = "Share Receiving", fill = "Race") 

ggsave("edu_share_transfer_income.png")

psid_edu %>%
  filter(Race_Head %in% c("White", "Black"),
         degree_type_final %in% c("2yr", "4yr"),
         !is.na(log_nonasset_income),
         Help_School_Indicator == 1,
         Help_School_Amount_Adj > 0) %>%
  mutate(income_quintile = ntile(log_nonasset_income, 4)) %>%
  group_by(Race_Head, income_quintile) %>%
  summarize(mean_amt = mean(Help_School_Amount_Adj, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(income_quintile), y = mean_amt, fill = Race_Head)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::dollar) +
  labs(title = "Mean Education Transfer Amount by Parental Income Quantile",
       subtitle = "Conditional on Receiving a Transfer (2010 $)",
       x = "Income Quantile (1 = Lowest)", y = "Mean Amount", fill = "Race") 

ggsave("edu_mean_transfer_income.png")

library(marginaleffects)

plot_predictions(pooled_probit_edu, 
                 condition = c("log_nonasset_income", "Race_Head")) +
  labs(title = "Predicted Probability of Education Transfer by Income",
       x = "Log Non-Asset Income", y = "Predicted P(Transfer)", color = "Race_Head") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("edu_transfer_income_prediction.png")

plot_predictions(pooled_probit_edu,
                       condition = c("Head_College", "Race_Head")) +
  labs(title = "By Parental Education",
       x = "Parental Education", y = "Predicted P(Transfer)") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA))


library(margins)

mfx <- margins(pooled_probit_edu)

summary(mfx) %>%
  as_tibble() %>%
  filter(factor %in% c("log_nonasset_income", "log_asset_income",
                       "degree_type_final4yr", "Head_CollegeNo College",
                       "Head_CollegeSome College", "Race_HeadWhite",
                       "enroll_erapre-2002", "enroll_era2009-2012")) %>%
  ggplot(aes(x = reorder(factor, AME), y = AME)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(title = "Average Marginal Effects: Probit for Transfer Receipt",
       x = NULL, y = "Average Marginal Effect on P(Transfer)") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA))

library(patchwork)

p1 <- plot_predictions(pooled_probit_edu, 
                       condition = c("log_nonasset_income", "Race_Head")) +
  labs(title = "By Labor Income",
       x = "Log Non-Asset Income", y = "Predicted P(Transfer)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

p2 <- plot_predictions(pooled_probit_edu, 
                       condition = c("log_asset_income", "Race_Head")) +
  labs(title = "By Asset Income",
       x = "Log Asset Income", y = "Predicted P(Transfer)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

p3 <- plot_predictions(pooled_probit_edu,
                       condition = c("Head_College", "Race_Head")) +
  labs(title = "By Parental Education",
       x = "Parental Education", y = "Predicted P(Transfer)") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA))

gridExtra::grid.arrange(p1, p2, p3, nrow = 1, 
             top = "Predicted Probability of Education Transfer")


ggsave("edu_transfer_predictions.png", width = 14, height = 5)