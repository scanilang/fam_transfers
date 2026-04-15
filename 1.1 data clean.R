library(dplyr)

#######################################################################################
# Read in raw data 
#######################################################################################
psid_fam<- read.csv("../data/famtransfer_famfile.csv") %>% rename(ER30001 = ID_1968)
psid_ind<- read.csv("../data/famtransfer_indfile.csv") %>% filter(ER30001 < 3000) %>% select(-Interview_Number)
cpi_data <- openxlsx::read.xlsx("../data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 

# function to take the mode
getmode <- function(v) {
  v = na.omit(v)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#######################################################################################
# Main PSID HH
#######################################################################################
psid_clean = psid_ind %>% 
  # limit to household head
  filter(Relationship_Head %in% c(1,10) & Sequence_Number == 1) %>% 
  left_join(psid_fam, multiple = 'all') %>% # merge individual and family file to get individual id
  mutate(Year = Survey_Year - 1) %>% 
  # take mode of reported Year Born
  group_by(ER30001, ER30002) %>% 
  mutate(Age = if_else(Age == 999, NA, Age),
         Race_Head = case_when(Race_Head == 1 ~ "White",
                               Race_Head== 2 ~ "Black",
                               TRUE ~ as.character(Race_Head)),
         Year_Born = if_else(Year_Born == 9999, NA, Year_Born),
         Year_Born = getmode(Year_Born),
         Year_Born = if_else(is.na(Age), Survey_Year - Age_Head, if_else(is.na(Year_Born), Survey_Year - Age,Year_Born))) %>% 
  ungroup() %>% 
  mutate(Age_recode = Year - Year_Born) %>% 
  # variable cleaning
  mutate(# Demographic Variables
         Age_Head = case_when(Age_Head == 99 & Survey_Year %in% 1985:1996 ~ NA,
                              Age_Head == 999 ~ NA,
                              Age_Head == 98 & Survey_Year %in% c(1995,1996) ~ NA, 
                              Age_Head == 0 & Survey_Year %in% c(1994,1995,1996) ~ NA,
                              TRUE ~ Age_Head),
         Age_Spouse = case_when(Age_Spouse ==99  & Survey_Year %in% 1985:1997 ~ NA,
                                Age_Spouse == 999 ~ NA,
                                TRUE ~ Age_Spouse),
         Marital_Status = case_when(Marital_Status == 1 ~ "Married",
                                    Marital_Status == 2 ~ "Single",
                                    TRUE ~ as.character(Marital_Status)),
         # College
         Head_College_Degree = case_when(Head_College_Degree == 1 ~ "College Degree",
                                         Head_Years_College == 5 ~ "Attended, No Degree",
                                         Head_College_Degree %in% 8:9 ~ NA),
         Head_Years_College = if_else(Head_Years_College == 9, NA, Head_Years_College),
         Head_College = case_when(Head_Highest_Degree %in% 2:8 ~ "College Degree",
                                  Head_Years_College > 0 & Head_Highest_Degree %in% 0:1 ~ "Some College",
                                  TRUE ~ "No College"),
         (across(c(School_Expenses, Other_School_Expenses), ~ case_when(. == c(999998,999999) ~ NA,
                                                                       TRUE ~ .))),
         Total_Education_Expenditure = if_else(Total_Education_Expenditure <0, 0, Total_Education_Expenditure),
         (across(c(Head_Last_Year_Attend_College, Spouse_Last_Year_Attend_College), ~ case_when(. %in% 97:99 ~ NA,
                                                                        TRUE ~ .))),
         Head_Enrolled = if_else((Head_Last_Year_Attend_College == Year |Head_Last_Year_Attend_College == 96 ) , 1, 0),
         Spouse_Enrolled = if_else((Spouse_Last_Year_Attend_College == Year |Spouse_Last_Year_Attend_College == 96 ), 1, 0),
         # Family Types
         Father_Education_Head = if_else(Father_Education_Head %in% c(9,99), NA, Father_Education_Head),
         Mother_Education_Head = if_else(Mother_Education_Head %in% c(9,99), NA, Mother_Education_Head),
         Father_Education_Spouse = if_else(Father_Education_Spouse %in% c(9,99), NA, Father_Education_Spouse),
         Mother_Education_Spouse = if_else(Mother_Education_Spouse %in% c(9,99), NA, Mother_Education_Spouse),
         parent_educ_head = pmax(Father_Education_Head,  Mother_Education_Head, na.rm = TRUE),
         parent_educ_spouse = pmax(Father_Education_Spouse, Mother_Education_Spouse, na.rm = TRUE),
         Parents_Poor_Head = if_else(Parents_Poor_Head %in% c(0, 8,9), NA, Parents_Poor_Head)) %>% 
  arrange(ER30001, ER30002, Survey_Year) %>%
  group_by(ER30001, ER30002) %>%
  tidyr::fill(Father_Education_Head, 
              Mother_Education_Head,
              Father_Education_Spouse, 
              Mother_Education_Spouse,
              Parents_Poor_Head,
              .direction = "downup") %>%
  ungroup() %>%
  mutate(# less than college, some college or nonacademic training, college or more
         head_side = case_when(
           # Primary: use Parents_Poor_Head when available (96% coverage)
           Parents_Poor_Head == 1 ~ "low",           # Poor
           Parents_Poor_Head == 3 ~ "mid",           # Average; "it varied"  
           Parents_Poor_Head == 5 ~ "high",          # Pretty well off
           
           # Fallback: use parental education when Parents_Poor_Head is NA
           parent_educ_head %in% 1:4 ~ "low",
           parent_educ_head %in% 5:6 ~ "mid",
           parent_educ_head >= 7     ~ "high",
           
           TRUE ~ NA_character_
         ),
         spouse_side = case_when(
           parent_educ_spouse %in% 1:4 ~ "low",
           parent_educ_spouse %in% 5:6       ~ "mid",
           parent_educ_spouse >= 7      ~ "high",
           TRUE ~ NA_character_),
         # family_type = case_when(
         #   # married households — both sides observable
         #   head_side == "high" & spouse_side == "high"                         ~ "both_high",
         #   (head_side == "high" & spouse_side %in% c("low", "mid")) |
         #     (head_side %in% c("low", "mid") & spouse_side == "high")          ~ "mixed",
         #   head_side == "mid"  & spouse_side == "mid"                          ~ "both_mid",
         #   head_side == "low"  & spouse_side == "low"                          ~ "both_low",
         #   head_side == "mid"  & spouse_side == "low"                          ~ "both_low",
         #   head_side == "low"  & spouse_side == "mid"                          ~ "both_low",
         #   
         #   # single households — spouse_side is NA, use head side only
         #   is.na(spouse_side) & head_side == "high"                            ~ "mixed",
         #   is.na(spouse_side) & head_side == "mid"                             ~ "both_mid",
         #   is.na(spouse_side) & head_side == "low"                             ~ "both_low",
         #   TRUE ~ NA_character_ ),
         family_type = case_when(
           head_side == "low"  ~ "low",
           head_side == "mid"  ~ "mid",
           head_side == "high" ~ "high",
           TRUE ~ NA_character_
         ),
         # Receive Support Variables
         (across(c(Received_Family_Reported_Head,Received_Family_Reported_Spouse,
                   Received_ChildSupport_Reported_Head, Received_ChildSupport_Reported_Spouse,
                   Received_Alimony_Reported_Head), 
                       ~case_when(. == 99999 & Survey_Year == 1993 ~ NA,
                                  . %in% c(99998,99999) & Survey_Year %in% 1994:2001 ~ NA,
                                  . %in% c(999998,999999) & Survey_Year %in% 2003:2021 ~ NA,
                                  TRUE ~ . ))),
         Received_Family_Head = if_else(is.na(Received_Family_Amount_Head),Received_Family_Reported_Head, Received_Family_Amount_Head),
         Received_Family_Spouse = if_else(is.na(Received_Family_Amount_Spouse),Received_Family_Reported_Spouse, Received_Family_Amount_Spouse),
         Received_ChildSupport_Head = if_else(is.na(Received_ChildSupport_Amount_Spouse),Received_ChildSupport_Reported_Head, Received_ChildSupport_Amount_Spouse),
         Received_ChildSupport_Spouse = if_else(is.na(Received_ChildSupport_Amount_Spouse),Received_ChildSupport_Reported_Spouse, Received_ChildSupport_Amount_Spouse),
         Received_Alimony_Head = if_else(is.na(Received_Alimony_Amount_Head), Received_Alimony_Reported_Head, Received_Alimony_Amount_Head),
         Received_ChildSupport_Alimony = Received_ChildSupport_Head + Received_ChildSupport_Spouse + Received_Alimony_Head,
         Received_Support_Amount_Head_Spouse = Received_Family_Head + Received_Family_Spouse, # this already does not include alimony and child support
         # 1988 Survey Synchronize 
         Supported1 = if_else(Survey_Year == 1988 & Provided_Support_ExPartner_88 > 0, "Ex/Current Partner", as.character(Supported1)),
         Supported2 = if_else(Survey_Year == 1988 & Provided_Support_Children_88 > 0, "Child", as.character(Supported2)),
         Supported3 = if_else(Survey_Year == 1988 & Provided_Support_Parents_88 > 0, "Parent", as.character(Supported3)),
         Supported4 = if_else(Survey_Year == 1988 & Provided_Support_Siblings_88 > 0, "Sibling", as.character(Supported4)),
         Supported5 = if_else(Survey_Year == 1988 & Provided_Support_OtherRelative_88 > 0, "Other relative", as.character(Supported5)),
         Provided_Family_Support_Total_88 = if_else(Survey_Year == 1988, Provided_Support_Total_88 - Provided_Support_NonRelative_88 - Provided_Support_Unknown_88, Provided_Support_Total_88), 
         # Provide Support Variables
         (across(c(Provided_Support_Amount_Head_Spouse, Provided_Alimony_Amount_Head_Spouse, Provided_ChildSupport_Amount_Head_Spouse), 
                 ~case_when(. == 99999 & Survey_Year <= 1992 ~ NA,
                            . == 9999999 & Survey_Year == 1993 ~ NA,
                            . %in% c(9999998,9999999) & Survey_Year %in% 1994:2021 ~ NA,
                            TRUE ~ . ))),
         (across(c(Supported1, Supported2, Supported3, Supported4, Supported5), 
                   ~case_when(. %in% 10:29 & Survey_Year != 1994 ~ "Ex/Current Partner",
                              . %in% 30:39 & Survey_Year != 1994 ~ "Child",
                              . %in% 40:49 & Survey_Year != 1994 ~ "Sibling",
                              . %in% 50:59 & Survey_Year != 1994 ~"Parent",
                              . %in% 60:65 & Survey_Year != 1994 ~ "(Great)grandchild",
                              . %in% 66:69 & Survey_Year != 1994 ~ "(Great)grandparent",
                              . %in% 70:97 & Survey_Year != 1994 ~ "Other relative",
                              . %in% 98 & Survey_Year != 1994 ~ "Nonrelative",
                              . %in% c(16,17,18,19) & Survey_Year == 1994  ~ "Ex or deceased spouse",
                              . %in% c(10,20,22,25,26) & Survey_Year == 1994  ~ "Partner",
                              . %in% 30:43 & Survey_Year == 1994 ~ "Child",
                              . %in% 45:54 & Survey_Year == 1994 ~ "Sibling",
                              . %in% 55:64 & Survey_Year == 1994 ~ "Parent",
                              . %in% c(65,66,71,72) & Survey_Year == 1994 ~ "(Great)grandchild",
                              . %in% c(67:70,73:76) & Survey_Year == 1994 ~ "(Great)grandparent",
                              . %in% 77:97 & Survey_Year == 1994 ~ "Other relative",
                              . %in% c(27,28) & Survey_Year == 1994 ~ "Nonrelative",
                              is.na(.) ~ "NA",
                             TRUE ~ as.character(.)))),
         Provided_Support_Nonrelative = if_else(if_any(starts_with("Supported"), ~ . == "Nonrelative"), 1, 0),
         Provided_Support_Child = if_else(if_any(starts_with("Supported"), ~ . == "Child"), 1, 0),
         Provided_Support_ExPartner = if_else(if_any(starts_with("Supported"), ~ . %in% c("Ex or deceased spouse","Ex/Current Partner")), 1, 0),
         Provided_Support_Unknown = if_else((if_any(starts_with("Supported"), ~ . %in% c("99","998","999"))) , 1, 0), 
         Provided_Support_Relative = if_else(if_any(starts_with("Supported"), ~ Provided_Support_Nonrelative != 1 & Provided_Support_Unknown != 1), 1, 0),
         Provided_ChildSupport_Amount_Head_Spouse = if_else(Provided_Support_Child == 0 &  Provided_Support_ExPartner == 0 , 0,Provided_ChildSupport_Amount_Head_Spouse),
         Provided_Family_Support_Total = case_when(Survey_Year == 1988 ~ Provided_Family_Support_Total_88,
                                                   Survey_Year != 1988 & Provided_Support_Amount_Head_Spouse == 0 ~ 0,
                                                   # Unidentified recipient 
                                                   (Survey_Year != 1988 & Provided_Support_Nonrelative == 1 & Provided_Support_Relative == 1 & 
                                                      Provided_Support_Amount_Head_Spouse- Provided_ChildSupport_Amount_Head_Spouse - Provided_Alimony_Amount_Head_Spouse> 0) ~ NA,
                                                   # Nonrelatives only
                                                   Survey_Year != 1988 & Provided_Support_Nonrelative ==1 & Provided_Support_Relative == 0 ~ 0,
                                                   Survey_Year != 1988 & Provided_Support_Unknown ==1 & Provided_Support_Relative == 0 ~ 0,
                                                   if_all(starts_with("Supported"), ~ . == "0") ~ 0,
                                                   TRUE ~ Provided_Support_Amount_Head_Spouse), 
         Provided_Family_Support_NetACS = Provided_Family_Support_Total- Provided_ChildSupport_Amount_Head_Spouse - Provided_Alimony_Amount_Head_Spouse,
         Provided_Family_Support_NetACS = if_else(Provided_Family_Support_NetACS < 0, 0 , Provided_Family_Support_NetACS), 
         # Labor Income not including Farm
         across(c(Labor_Income_LessFB_Reported_Spouse), ~if_else(. %in% c(9999998,9999999), NA, .)),
         Labor_Income_Spouse = if_else(is.na(Labor_Income_LessFB_Amount_Spouse), Labor_Income_LessFB_Reported_Spouse, Labor_Income_LessFB_Amount_Spouse),
         Labor_Income_Plus_Business_Spouse = case_when(Survey_Year < 1993 ~ Labor_Income_Amount_Spouse, # no farm labor income for spouse before 1993
                                                       Survey_Year == 1993 ~ Labor_Income_Amount_Spouse - Farm_Labor_Income_Spouse, 
                                                       Survey_Year >= 1994 ~ Labor_Income_Spouse +  Business_Income_Spouse), 
         Labor_Income_Plus_Business_Head = if_else(Survey_Year < 1994,Labor_Income_Amount_Head - Farm_Labor_Income_Head, 
                                                   Labor_Income_LessFB_Amount_Head + Business_Income_Head),
         Total_Labor_Income_Plus_Business = Labor_Income_Plus_Business_Head + Labor_Income_Plus_Business_Spouse,
         across(c(Unemployment_Comp_Head_Reported,Unemployment_Comp_Spouse_Reported,
                  Workers_Comp_Head_Reported, Workers_Comp_Spouse_Reported), ~case_when(. == 99999 & Survey_Year %in% 1993:2001 ~ NA,
                                                                                        . == 99998 & Survey_Year %in% 1994:2001 ~ NA,
                                                                                        . >= 999998 & Survey_Year %in% 2003:2021 ~ NA,
                                                                                        TRUE ~ .)),
         Unemployment_Comp_Head = if_else(is.na(Unemployment_Comp_Head),Unemployment_Comp_Head_Reported, Unemployment_Comp_Head),
         Unemployment_Comp_Spouse= if_else(is.na(Unemployment_Comp_Spouse),Unemployment_Comp_Spouse_Reported, Unemployment_Comp_Spouse),
         Workers_Comp_Head = if_else(is.na(Workers_Comp_Head),Workers_Comp_Head_Reported, Workers_Comp_Head),
         Workers_Comp_Spouse = if_else(is.na(Workers_Comp_Spouse),Workers_Comp_Spouse_Reported, Workers_Comp_Spouse),
         Total_UnempWorkers_Comp = Unemployment_Comp_Head + Unemployment_Comp_Spouse + Workers_Comp_Head + Workers_Comp_Spouse,
         # Asset Income not including farm 
         Asset_Income_Dividends_Reported_Head_raw = Asset_Income_Dividends_Reported_Head,
         (across(c(Asset_Income_Rent_Reported_Head, Asset_Income_TrustFund_Reported_Head, Asset_Income_Dividends_Reported_Head, Asset_Income_Interest_Reported_Head), 
                 ~case_when(. %in% c(-99999,999999) & Survey_Year == 1993 ~ NA,
                            . %in% c(-99999,999998, 999999) & Survey_Year %in% 1994:2021 ~ NA,
                            TRUE ~ . ))),
         (across(c(Asset_Income_TrustFund_Reported_Spouse, Asset_Income_Interest_Reported_Spouse, Asset_Income_Dividends_Reported_Spouse), 
                 ~case_when(. %in% c(99999) & Survey_Year == 1993 ~ NA,
                            . %in% c(-99999,999998, 999999) & Survey_Year %in% 1994:2021 ~ NA,
                            TRUE ~ . ))),
         Asset_Income_Rent_Head = if_else(is.na(Asset_Income_Rent_Amount_Head),Asset_Income_Rent_Reported_Head, Asset_Income_Rent_Amount_Head),
         Asset_Income_TrustFund_Head = if_else(is.na(Asset_Income_TrustFund_Amount_Head),Asset_Income_TrustFund_Reported_Head, Asset_Income_TrustFund_Amount_Head),
         Asset_Income_Interest_Head = if_else(is.na(Asset_Income_Interest_Amount_Head),Asset_Income_Interest_Reported_Head, Asset_Income_Interest_Amount_Head),
         Asset_Income_Dividends_Head = if_else(is.na(Asset_Income_Dividends_Amount_Head),Asset_Income_Dividends_Reported_Head, Asset_Income_Dividends_Amount_Head),
         Asset_Income_TrustFund_Spouse= if_else(is.na(Asset_Income_TrustFund_Amount_Spouse),Asset_Income_TrustFund_Reported_Spouse, Asset_Income_TrustFund_Amount_Spouse),
         Asset_Income_Interest_Spouse = if_else(is.na(Asset_Income_Interest_Amount_Spouse),Asset_Income_Interest_Reported_Spouse, Asset_Income_Interest_Amount_Spouse),
         Asset_Income_Dividends_Spouse = if_else(is.na(Asset_Income_Dividends_Amount_Spouse),Asset_Income_Dividends_Reported_Spouse, Asset_Income_Dividends_Amount_Spouse),
         Asset_Income_Head = Asset_Income_Dividends_Head + Asset_Income_Interest_Head + Asset_Income_TrustFund_Head + Asset_Income_Business_Head + Asset_Income_Rent_Head,
         Asset_Income_Spouse = Asset_Income_Dividends_Spouse + Asset_Income_Interest_Spouse + Asset_Income_TrustFund_Spouse + Asset_Income_Business_Spouse,
         Total_Asset_Income = if_else(Survey_Year < 1993, Asset_Income_Business_Head_Spouse +Asset_Income_All_Head + Asset_Income_All_Spouse,
                                      Asset_Income_Head + Asset_Income_Spouse),
         # Total_Asset_Income = Total_Taxable_Income_Head_Spouse - Labor_Income
         # Public Transfers (not including social security)
         (across(c(ADC_AFDC_Reported_Head, SSI_Reported_Head, Other_Welfare_Reported_Head), 
                 ~case_when(. == 99999 & Survey_Year == 1993 ~ NA,
                            . %in% c(999998,999999) & Survey_Year %in% 1994:2021 ~ NA,
                            TRUE ~ . ))),
         ADC_AFDC_Head = if_else(is.na(ADC_AFDC_Amount_Head), ADC_AFDC_Reported_Head, ADC_AFDC_Amount_Head),
         SSI_Head = if_else(is.na(SSI_Amount_Head), SSI_Reported_Head, SSI_Amount_Head),
         Other_Welfare_Head = if_else(is.na(Other_Welfare_Amount_Head), Other_Welfare_Reported_Head, Other_Welfare_Amount_Head),
         Public_Transfers_Head = ADC_AFDC_Head + SSI_Head + Other_Welfare_Head,
         (across(c(ADC_AFDC_Reported_Spouse, SSI_Reported_Spouse, Other_Welfare_Reported_Spouse), 
                 ~case_when(. == 99999 & Survey_Year == 1993 ~ NA,
                            . %in% c(99998,99999) & Survey_Year %in% 1994:2001 ~ NA,
                            . %in% c(999998,999999) & Survey_Year %in% 2003:2021 ~ NA,
                            TRUE ~ . ))),
         ADC_AFDC_Spouse = if_else(is.na(ADC_AFDC_Amount_Spouse), ADC_AFDC_Reported_Spouse, ADC_AFDC_Amount_Spouse),
         SSI_Spouse = if_else(is.na(SSI_Amount_Spouse), SSI_Reported_Spouse, SSI_Amount_Spouse),
         Other_Welfare_Spouse = if_else(is.na(Other_Welfare_Amount_Spouse), Other_Welfare_Reported_Spouse, Other_Welfare_Amount_Spouse),
         Public_Transfers_Spouse = ADC_AFDC_Spouse + SSI_Spouse + Other_Welfare_Spouse,
         Total_Public_Transfers = case_when(Survey_Year == 1985 ~ ADC_AFDC_Head_Spouse + SSI_Head_Spouse + Other_Welfare_Head_Spouse,
                                            Survey_Year %in% c(1986:1993, 2005:2021) ~ Public_Transfers_Head + Public_Transfers_Spouse,
                                            Survey_Year %in% 1994:2004 & Age_recode < 62 &(Age_Spouse < 62| Age_Spouse == 0) ~ Public_Transfers_Head + Public_Transfers_Spouse,
                                            Survey_Year %in% 1994:2004 & Age_recode >= 62 | (Age_Spouse >= 62 | Age_Spouse == 0)~ Public_Transfers_Head + Public_Transfers_Spouse),
         # Retirement Income (including social security but excluding welfare)
         Social_Security_Total_Family = if_else(Social_Security_Total_Family == 999999 & Survey_Year %in% 1994:1995, NA,Social_Security_Total_Family ),
         (across(c(Retire_Pension_Reported_Head, Annuity_IRA_Reported_Head, Other_Retirement_Reported_Head), 
                 ~case_when(. == 99999 & Survey_Year == 1993 ~ NA,
                            . %in% c(999998,999999) & Survey_Year %in% 2003:2021 ~ NA,
                            TRUE ~ . ))),
         Retire_Pension_Head = if_else(is.na(Retire_Pension_Amount_Head), Retire_Pension_Reported_Head, Retire_Pension_Amount_Head),
         Annuity_IRA_Head = if_else(is.na(Annuity_IRA_Amount_Head), Annuity_IRA_Reported_Head, Annuity_IRA_Amount_Head),
         Other_Retirement_Head = if_else(is.na(Other_Retirement_Amount_Head), Other_Retirement_Reported_Head, Other_Retirement_Amount_Head),
         (across(c(Other_Retirement_All_Reported_Spouse, VA_Pension_Reported_Head), 
                 ~case_when(. == 99999 & Survey_Year == 1993 ~ NA,
                            . %in% c(99998,99999) & Survey_Year %in% 1994:2001 ~ NA,
                            . %in% c(999998,999999) & Survey_Year %in% 2003:2021 ~ NA,
                            TRUE ~ . ))),
         Other_Retirement_All_Spouse = if_else(is.na(Other_Retirement_All_Amount_Spouse), Other_Retirement_All_Reported_Spouse, Other_Retirement_All_Amount_Spouse),
         VA_Pension_Head = if_else(is.na(VA_Pension_Amount_Head), VA_Pension_Reported_Head, VA_Pension_Amount_Head),
         Total_Retirement_Income = case_when(Survey_Year < 1993 ~ Other_Retirement_All_Head + Other_Retirement_All_Spouse + VA_Pension_Head,
                                             Survey_Year %in% 1993:2011 ~ Retire_Pension_Head + Annuity_IRA_Head + Other_Retirement_Head +VA_Pension_Head + Other_Retirement_All_Spouse,
                                             Survey_Year %in% 2013:2021 ~ Retire_Pension_Head + Annuity_IRA_Head + Other_Retirement_Head +VA_Pension_Head + Other_Retirement_Spouse + Annuities_Spouse + IRAs_Spouse +Retirement_Pension_Spouse ),
         Total_Retirement_Income = case_when(Survey_Year == 1985 ~ Total_Retirement_Income + Social_Security_Head_Spouse,
                                             Survey_Year %in% c(1986:1993, 2005:2021) ~ Total_Retirement_Income + Social_Security_Head + Social_Security_Spouse,
                                             Survey_Year %in% 1994:2004 & Age_recode < 62 &(Age_Spouse < 62| Age_Spouse == 0) ~ Total_Retirement_Income,
                                             Survey_Year %in% 1994:2004 & Age_recode >= 62 | (Age_Spouse >= 62 | Age_Spouse == 0)~ Total_Retirement_Income+ Social_Security_Total_Family),
         # Summary
         Labor_UIWC_Income = Total_Labor_Income_Plus_Business + Total_UnempWorkers_Comp, 
         Family_Transfers_Net = Received_Support_Amount_Head_Spouse - Provided_Family_Support_NetACS,
         Total_Income_Head_Spouse = (Total_Labor_Income_Plus_Business + Total_Asset_Income +  Total_Public_Transfers + Total_UnempWorkers_Comp + Received_ChildSupport_Alimony +
                                       Total_Retirement_Income- (Provided_ChildSupport_Amount_Head_Spouse - Provided_Alimony_Amount_Head_Spouse)),
         ChildSupport_Alimony_Adjustment = Received_ChildSupport_Alimony- (Provided_ChildSupport_Amount_Head_Spouse + Provided_Alimony_Amount_Head_Spouse),
         Total_NonAsset_Income = Total_UnempWorkers_Comp + Total_Labor_Income_Plus_Business + Total_Public_Transfers+ Total_Retirement_Income + ChildSupport_Alimony_Adjustment) %>% 
  # Attempt to figure out Youngest Child in or out of FU
  group_by(ER30001, ER30002) %>% 
  mutate(Most_Recent_Child_Born = if_else(Child_Born_Head_Spouse %in% 1:3| Child_Born_Head_Only %in% 1:3 | Child_Born_Spouse_Only %in% 1:3, Year, NA_integer_),
         Most_Recent_Child_Born = if_else(Child_Born_Head_Spouse_2Years %in% 1:3| Child_Born_Head_Only_2Years %in% 1:3 | Child_Born_Spouse_Only_2Years %in% 1:3, Year-1, Most_Recent_Child_Born)) %>% 
  tidyr::fill(Most_Recent_Child_Born, .direction = "down" ) %>% 
  mutate(Age_Youngest = Year - Most_Recent_Child_Born, 
         Have_Dependents = if_else(Age_Youngest <= 21 | Number_Dependents_Head_Spouse > 0, 1, 0)) %>% 
  ungroup() %>% 
  # CPI adjustment
  left_join(cpi_data %>% select(year, ratio_2010), by = c("Year" = "year")) %>% 
  mutate(across(c(ChildSupport_Alimony_Adjustment,Received_Support_Amount_Head_Spouse, Provided_Family_Support_NetACS, Family_Transfers_Net,
                  Total_Retirement_Income,Total_Asset_Income, Total_Income_Head_Spouse,Labor_UIWC_Income,
                  Total_Labor_Income_Plus_Business, Total_Public_Transfers , Unemployment_Comp_Head, Workers_Comp_Head, 
                  Unemployment_Comp_Spouse, Workers_Comp_Spouse, Total_NonAsset_Income,
                  Total_Family_Income, Total_Taxable_Income_Head_Spouse, Total_Transfer_Head_Spouse, Total_Taxable_Income_Other,
                  Total_Transfers_Other, Total_Social_Security_Other, School_Expenses, Other_School_Expenses, Total_Education_Expenditure), ~ . *ratio_2010)) %>% 
  arrange(ER30001, ER30002, Survey_Year) %>%
  group_by(ER30001, ER30002) %>%
  mutate(received_transfer_past = cumsum(Received_Support_Amount_Head_Spouse > 0) > 0 & row_number() > 1,
         provided_transfer_past = cumsum(Provided_Family_Support_NetACS > 0) > 0 & row_number() > 1) %>% 
  ungroup() %>% 
  # Variables for transfer functions
  mutate(Provided_Support_Indicator = if_else(Provided_Family_Support_NetACS > 0, 1, 0),
         Received_Support_Indicator = if_else(Received_Support_Amount_Head_Spouse > 0, 1, 0),
         log_asset_income = log(Total_Asset_Income + 1),
         log_nonasset_income = log(Total_NonAsset_Income + 1),
         log_transfer_in = log(Received_Support_Amount_Head_Spouse + 1),
         log_transfer_out = log(Provided_Family_Support_NetACS + 1),
         Age2 = Age^2,
         # Year_bins = paste(1983 + ((Year-1983) %/% 10) * 10, 
         #                 1992 + ((Year-1983) %/% 10) * 10, 
         #                 sep = "-"),
         # Year_bins = forcats::fct_relevel(Year_bins, c("2013-2022", "2003-2012", "1993-2002", "1983-1992")),
         Year_bins = case_when(
           Year %in% 1984:1993 ~ "1984-1993",   # includes 90-91 recession
           Year %in% 1994:2007 ~ "1994-2007",   # long expansion
           Year %in% 2008:2012 ~ "2008-2012",   # financial crisis and aftermath
           Year %in% 2013:2021 ~ "2013-2021",   # recovery and COVID
           TRUE ~ NA_character_
         ),
       #   Birth_Cohort = paste((Year_Born %/% 10) * 10, 
       #                      (Year_Born %/% 10) * 10 + 9, 
       #                      sep = "-"),
       #   # regrouping first and last group
       # Birth_Cohort = case_when(Birth_Cohort == "2000-2009"~ "1990-1999",
       #                          Birth_Cohort == "1910-1919"~ "1920-1929", 
       #                          Birth_Cohort == "1900-1909"~ "1920-1929", 
       #                          TRUE ~ Birth_Cohort),
       # Birth_Cohort = forcats::fct_relevel(Birth_Cohort, c("1970-1979", "1980-1989", "1990-1999", "1960-1969", "1950-1959","1940-1949","1930-1939", "1920-1929")),
         Birth_Cohort = case_when(
           Year_Born < 1940              ~ "pre-1940",    # collapse thin early cohorts
           Year_Born %in% 1940:1954      ~ "1940-1954",
           Year_Born %in% 1955:1969      ~ "1955-1969",
           Year_Born %in% 1970:1984      ~ "1970-1984",
           Year_Born >= 1985             ~ "1985-plus",
         TRUE ~ NA_character_
       ),
       log_labor_uiwc = log(Labor_UIWC_Income + 1))
  

write.csv(psid_clean, "../data/psid_clean.csv")

# 
# View(psid_clean %>% select(ER30001, ER30002, Survey_Year, Marital_Status, Youngest_in_FU, Num_Children_FU, Number_Dependent_Outside_FU, Provided_ChildSupport_Amount_Head_Spouse, Provided_Alimony_Amount_Head_Spouse,
#                            Received_ChildSupport_Spouse, Received_ChildSupport_Head , Received_Alimony_Head, Age_Youngest,Have_Dependents, rec_al, rec_cs, prov_al, prov_cs) %>% 
#        filter((prov_cs > 0 | prov_al > 0 | rec_cs > 0 | rec_al > 0)))
# 
# View(psid_clean %>% select(ER30001, ER30002, Survey_Year, Asset_Income_Dividends_Reported_Head_raw, Asset_Income_Dividends_Head, Asset_Income_Interest_Head , Asset_Income_TrustFund_Head , Asset_Income_Business_Head ,
#                            Asset_Income_Rent_Head,Asset_Income_Dividends_Spouse ,Asset_Income_Interest_Spouse , Asset_Income_TrustFund_Spouse , Labor_Income_Plus_Business_Head,
#                            Asset_Income_Business_Spouse,Asset_Income_Business_Head_Spouse ,Asset_Income_All_Head , Asset_Income_All_Spouse, Total_Asset_Income) %>% filter(is.na(Total_Asset_Income)))
# ###################################################
# Testing 
###################################################

# Receive_Family_Reported_Head:
#topcoded to 999997 from 2003 to 2021
#topcoded to 99997 from 1994 to 2001
#topcoded 99,998 in 1993

# Provide_Support_Head:
#topcoded to 99998 from 1985-1992
#topcoded to 9999998 in 1993
#topcoded to 9999997 in 1994-2021

# How much is filtered out when excluding all that answered non_relative
# 318 gave to both a nonrelative and listed relative (non refused answer)
# test= psid_clean %>% 
#   mutate(Provided_Support_Relative = if_else(if_any(starts_with("Supported"), ~ . != "Nonrelative" & .!= "0" & .!= "99" & .!= "999" & .!= "998"), 1, 0)) %>% 
#   filter(Provided_Support_Nonrelative == 1 & Provided_Support_Relative ==1 ) %>% 
#   filter(Provided_Support_NetACS > 0) %>% 
#   select(Family_ID, ER30001,Year, Survey_Year,Age_Head, Race_Head, Marital_Status, Num_Children_FU,Provided_Support_Amount_Head_Spouse, Provided_Alimony_Amount_Head_Spouse,
#          Provided_ChildSupport_Amount_Head_Spouse,Provided_Support_NetACS, Supported1, Supported2, Supported3, Supported4, Supported5, Number_Supported)
