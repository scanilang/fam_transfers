library(dplyr)

###################################################
# Read in clean data 
###################################################
psid_clean = read.csv("../data/psid_clean.csv")

###################################################
# Aggregate to 2 years
###################################################

nonaggregate_data <- psid_clean %>% 
  left_join(cpi_data %>% select(year, ratio_2010), by = c("Year" = "year")) %>% 
  select(ER30001, ER30002, Relationship_Head,Family_ID,Survey_Year, Year,Year_Born, Year_Born_bin, Sex, Race,Hispanic, Youngest_in_FU, ratio_2010, Age,Completed_Education,  Num_Children_FU, Family_Unit_Size, 
         Sex, Marital_Status, Received_Family_Indicator, Provide_Support_Head, Receive_Head_Spouse_Net, Support_Amount_Net,  
         Labor_Income_Household, Labor_Income_Head, Labor_Income_Spouse, Total_Family_Income,Total_Labor_Income, Unemployment_Comp_Head,Unemployment_Comp_Spouse, Workers_Comp_Head,Workers_Comp_Spouse,
         Supported1, Supported2, Supported3, Supported4, Supported5) %>% 
  mutate(Total_Family_Income = Total_Family_Income *ratio_2010, 
         Receive_Head_Spouse_Net = Receive_Head_Spouse_Net * ratio_2010,
         Support_Amount_Net = Support_Amount_Net * ratio_2010,
         Unemployment_Comp_Head  = Unemployment_Comp_Head * ratio_2010,
         Workers_Comp_Head  = Workers_Comp_Head * ratio_2010,
         Unemployment_Comp_Spouse  = Unemployment_Comp_Spouse * ratio_2010,
         Workers_Comp_Spouse  = Workers_Comp_Spouse * ratio_2010,
         Total_Labor_Income = Total_Labor_Income *ratio_2010,
         Labor_Income_Household =Labor_Income_Household *ratio_2010,
         Labor_Income_Head =Labor_Income_Head *ratio_2010,
         Labor_Income_Spouse =Labor_Income_Spouse *ratio_2010) %>% 
  select(-ratio_2010) %>% 
  mutate(log_fam_inc = log(Total_Family_Income)) %>% 
  filter(Age>= 25,
         Race %in% c("Black", "White"),
         Hispanic %in% c(NA, "Inap")) %>% 
  filter(Marital_Status %in% c("Married", "Single")) 

# aggregate 4 years
data_4yr_aggregate <- nonaggregate_data %>% 
  select(ER30001, ER30002, Survey_Year, Year, Year_Born,Year_Born_bin,  Age,Sex, Race, Marital_Status,Family_Unit_Size,Num_Children_FU, Youngest_in_FU, Completed_Education,
         Support_Amount_Net, Receive_Head_Spouse_Net, Total_Family_Income, Total_Labor_Income,Labor_Income_Household, Labor_Income_Head, Labor_Income_Spouse, Unemployment_Comp_Head,Unemployment_Comp_Spouse, Workers_Comp_Head,Workers_Comp_Spouse,
         Supported1, Supported2, Supported3, Supported4, Supported5) %>% 
  full_join(data_final %>% # double biennial year
            filter(Survey_Year >= 1999) %>% 
            mutate(Year = Year -1,Age = Age -1) %>% 
            select(ER30001, ER30002, Survey_Year, Year, Year_Born,Year_Born_bin, Age,Sex, Race, Marital_Status,Completed_Education,Family_Unit_Size,Num_Children_FU,
                   Youngest_in_FU, Support_Amount_Net, Receive_Head_Spouse_Net, Total_Family_Income,Total_Labor_Income,Labor_Income_Household, Labor_Income_Head, Labor_Income_Spouse, Unemployment_Comp_Head,Unemployment_Comp_Spouse, Workers_Comp_Head,Workers_Comp_Spouse,
                   Supported1, Supported2, Supported3, Supported4, Supported5)) %>% 
  arrange(ER30001, ER30002, Survey_Year, Year) %>% 
  # some manual edits to age
  mutate(Age = case_when(ER30001 == 1822 & ER30002 == 1 & Age == 79 ~ 81,
                         ER30001 == 1822 & ER30002 == 1 & Age == 80 ~ 82,
                         ER30001 == 1822 & ER30002 == 1 & Age == 81 ~ 83,
                         ER30001 == 1822 & ER30002 == 1 & Age == 82 ~ 84,
                         TRUE ~ Age)) %>% 
  mutate(age_cat = case_when(Age %in% 25:28 ~ 1,
                             Age %in% 29:32 ~ 2,
                             Age %in% 33:36 ~ 3,
                             Age %in% 37:40 ~ 4,
                             Age %in% 41:44 ~ 5,
                             Age %in% 45:48 ~ 6,
                             Age %in% 49:52 ~ 7,
                             Age %in% 53:56 ~ 8,
                             Age %in% 57:60 ~ 9,
                             Age %in% 61:64 ~ 10,
                             Age %in% 65:68 ~ 11,
                             Age %in% 69:72 ~ 12,
                             Age %in% 73:76 ~ 13,
                             Age %in% 77:80 ~ 14,
                             Age %in% 81:84 ~ 15,
                             Age %in% 85:88 ~ 16,
                             Age %in% 89:92 ~ 17)) %>% 
  group_by(ER30001, ER30002, age_cat) %>% 
    mutate(n = n(),
           sum_transfer_out = sum(Support_Amount_Net),
           sum_transfer_in = sum(Receive_Head_Spouse_Net),
           sum_labor_income_hh = sum(Labor_Income_Household),
           sum_labor_income_head = sum(Labor_Income_Head),
           sum_fam_income = sum(Total_Family_Income),
           sum_unemp_insurance = sum(Unemployment_Comp_Head + Unemployment_Comp_Spouse),
           sum_workers_comp = sum(Workers_Comp_Head + Workers_Comp_Spouse),
           sum_labor_income = sum(Total_Labor_Income),
           years_married = sum(Marital_Status == "Married"),
           avg_FU_size= mean(Family_Unit_Size),
           Completed_Education_max = max(Completed_Education),
           Year_first = min(Year),
           sum_labor_income_comp_ui = sum_labor_income_hh + sum_unemp_insurance + sum_workers_comp) %>% 
    ungroup() %>% 
   filter(n == 4) %>% 
  mutate(Marital_Status = if_else(years_married >= 3, "Married","Single"),
         age_cat = if_else(age_cat >=16, 15, age_cat)) %>% # for sample size, include any consecutive 4 non overlapping 81-92
 distinct(ER30001, ER30002, age_cat,Year_Born,Year_Born_bin,  Year_first,Race, Sex,Marital_Status, years_married,avg_FU_size,
          sum_fam_income, sum_labor_income,sum_labor_income_hh, sum_labor_income_head,sum_transfer_out,sum_transfer_in, sum_unemp_insurance, sum_workers_comp, 
          sum_labor_income_comp_ui, sum_asset_income, sum_asset_labor_income, sum_wealth) %>% 
 # filter out outliers
filter(!(ER30001 == 2225 & ER30002 == 4 & age_cat == 15 ),
       !(ER30001 == 1824 & ER30002 == 181 & age_cat == 11 )) %>% 
  mutate(Year_first_bin = case_when(Year_first %in% 1984:1987 ~ '1984-1990',
                                    Year_first %in% 1988:1997 ~ '1988-2000',
                                    Year_first %in% 1998:2007 ~ '1998-2010',
                                    Year_first %in% 2008:2017 ~ '2008-2020')) %>% 
  group_by(ER30001, ER30002) %>%
  mutate(received_transfer_past = cumsum(sum_transfer_in > 0) > 0 & row_number() > 1,
         provided_transfer_past = cumsum(sum_transfer_out > 0) > 0 & row_number() > 1) %>% 
  ungroup()

# write output
write.csv(data_4yr_aggregate, "PSID data/data_4yr_aggregate.csv")

View(data_4yr_aggregate %>% group_by(Race, age_cat)  %>%  count())

# Year, Race,Hispanic, ratio_2010, Age,Youngest_in_FU, Num_Children_FU, Family_Unit_Size, 
# Sex, Marital_Status, Received_Family_Indicator, Provide_Support_Head, Total_Family_Income,
# Wealth_W_Equity, Wealth_WO_Equity, Receive_Head_Spouse_Net, Support_Amount_Outside_FU, Support_Amount_Net,  Inheritance_Val_1, Inheritance_Val_2, Inheritance_Val_3,Inheritance_Indicator,
# Supported1, Supported2, Supported3, Supported4, Supported5)

View(psid_clean%>% select(ER30001, ER30002, Year,Total_Labor_Income,Wages_Head, Wages_Spouse, Unemployment_Comp_Head,Unemployment_Comp_Spouse, 
                          Workers_Comp_Head,Workers_Comp_Spouse,Employment_Status_Head,Total_Labor_Income_Head_FB))


         