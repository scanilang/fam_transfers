library(dplyr)

#################################################################################################################################
# Read in raw or cleaned data 
#################################################################################################################################
psid_ind<- read.csv("../data/famtransfer_indfile.csv") %>% filter(ER30001 < 3000) %>% select(-Interview_Number)
cpi_data <- openxlsx::read.xlsx("../data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 

load(file= "/Users/scanilang/Documents/econ/psidR/data/RT13PARCHD.rda")
psid_fam_roster <- get("RT13PARCHD")

# cleaned psid data
psid_clean <- read.csv("../data/psid_clean.csv")

# function to take the mode
getmode <- function(v) {
  v = na.omit(v)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#################################################################################################################################
# Enrollment Start and End Dates
#################################################################################################################################
college_enrollment <- psid_ind %>%
  select(ER30001, ER30002, Family_ID, Survey_Year, 
         Relationship_Head, Age, Completed_Education) %>% 
  mutate(Year = Survey_Year - 1,
         Age = case_when(Age %in% c(99, 999) ~ NA_real_,
                         Age == 0             ~ NA_real_,
                         TRUE                 ~ as.numeric(Age)),
         Completed_Education = if_else(
           Completed_Education %in% c(98, 99, 0), NA_real_,
           as.numeric(Completed_Education))) %>%
  arrange(ER30001, ER30002, Survey_Year) %>%
  group_by(ER30001, ER30002) %>%
  tidyr::fill(Completed_Education, .direction = "down" ) %>% 
  mutate(
    Completed_Education = cummax(tidyr::replace_na(Completed_Education, 0)),
    Completed_Education = if_else(Completed_Education == 0, NA_real_, Completed_Education),
    lag_edu = lag(Completed_Education),
    max_education = max(Completed_Education, na.rm = TRUE),
    max_education = if_else(is.infinite(max_education), NA_real_, max_education),
    
    # censoring flags
    last_obs_year = max(Year, na.rm = TRUE),
    edu_at_last_obs = if_else(Year == last_obs_year, Completed_Education, NA_real_),
    lag_edu_at_last = if_else(Year == last_obs_year, lag_edu, NA_real_),
    
    # enrollment start
    first_obs_already_enrolled = is.na(lag_edu) & Completed_Education > 12,
    entered_college = !is.na(lag_edu) & lag_edu <= 12 & Completed_Education > 12,
    entered_college_btw_waves = !is.na(lag_edu) & lag_edu <= 12 & Completed_Education > 13, # 2 year jump bc biannual survey
    
    year_enrolled_start = case_when(
      entered_college & !entered_college_btw_waves  ~ Year,
      (first_obs_already_enrolled | entered_college_btw_waves) ~ Year - (Completed_Education - 12),
      TRUE                       ~ NA_real_
    ),
    
    # graduation thresholds
    reached_2yr = Completed_Education >= 14,
    reached_4yr = Completed_Education >= 16,
    year_at_2yr = if_else(reached_2yr, Year, NA_real_),
    year_at_4yr = if_else(reached_4yr, Year, NA_real_)
  ) %>% 
  summarize(
    
    max_education          = first(max_education),
    enroll_start_year      = min(year_enrolled_start, na.rm = TRUE),
    last_obs_year     = first(last_obs_year),
    first_obs_year         = min(Year, na.rm = TRUE),
    year_at_max_edu        = min(Year[Completed_Education == max_education], na.rm = TRUE),
  
    # education was still rising when they left the survey
    edu_still_rising  = any(edu_at_last_obs > lag_edu_at_last & 
                              edu_at_last_obs > 12 & edu_at_last_obs <= 16, 
                            na.rm = TRUE),
    # age when they reached max education
    age_at_max_edu    = min(Age[Completed_Education == max_education], na.rm = TRUE),
    
    grad_year_2yr          = min(year_at_2yr, na.rm = TRUE),
    grad_year_4yr          = min(year_at_4yr, na.rm = TRUE),
    enrollment_imputed     = any(first_obs_already_enrolled, na.rm = TRUE),
    enrollment_unobserved = enrollment_imputed & !is.na(year_at_max_edu) & (first_obs_year == year_at_max_edu),
   
     .groups = "drop"
  ) %>%
  mutate(
    across(c(enroll_start_year, grad_year_2yr, grad_year_4yr, year_at_max_edu, age_at_max_edu),
           ~ if_else(is.infinite(.), NA_real_, .)),
    degree_type = case_when(
      edu_still_rising & max_education %in% 13:15 ~ NA_character_,  # might still be in school
      max_education %in% 13:15 ~ "2yr",
      max_education >= 16      ~ "4yr",
      TRUE                     ~ NA_character_
    ),
    enroll_end_year = case_when(
      enrollment_unobserved ~ NA_real_,
      degree_type == "4yr"  ~ grad_year_4yr,
      degree_type == "2yr"  ~ grad_year_2yr,
      TRUE                  ~ NA_real_
    )
  ) %>% 
  select(ER30001, ER30002, enroll_start_year, enroll_end_year, degree_type, age_at_max_edu, enrollment_imputed, enrollment_unobserved, edu_still_rising)

#################################################################################################################################
# Parent Income at Enrollment Date
#################################################################################################################################
psid_income <- psid_clean %>%
  select(
    ER30001, ER30002, Family_ID, Survey_Year,Year,
    Race_Head, Age_Head, Marital_Status,
    Head_College, Family_Unit_Size,
    log_nonasset_income, log_asset_income,
    Total_NonAsset_Income, Total_Asset_Income,
    family_type,
    Father_Education_Head, Mother_Education_Head
  )

#################################################################################################################################
# PSID Family Roster and Transfers
#################################################################################################################################

psid_rt13_clean <- psid_fam_roster %>%
  transmute(
    # ── Keys ──
    RT_Interview_Number      = RT13V61,
 #   Head_Sequence_Number     = RT13V62,
    Head_1968_ID             = RT13V63,
    Head_Person_Number       = RT13V64,
    Head_Age                 = RT13V65,
  #  Wife_1968_ID             = RT13V67,
   # Wife_Person_Number       = RT13V68,
    Record_Type              = RT13V70,
    Record_Number            = RT13V72,
    # ── Child identifiers ──
    Child_1968_ID            = RT13V74,
    Child_Person_Number      = RT13V75,
    Child_In_FU              = RT13V76,
    Relationship_to_Head     = RT13V77,
 #   Relationship_to_Wife     = RT13V78,
    # ── Child demographics ──
    Child_Birth_Year         = RT13V81,
    Child_Age                = RT13V83,
    Child_Education          = RT13V90,
    Child_Marital_Status     = RT13V91,
    Child_Has_Partner        = RT13V92,
    Child_Num_Kids           = RT13V93,
    # ── Money transfers (2012) ──
    Wtr_Parent_Gave_Money    = RT13V124,
    Amt_Parent_Gave_Money    = RT13V125,
    Amt_Parent_Gave_Per      = RT13V126,
    Wtr_Child_Gave_Money     = RT13V127,
    Amt_Child_Gave_Money     = RT13V128,
    Amt_Child_Gave_Per       = RT13V129,
    # ── Cumulative transfers (since 18) ──
    Amt_School_Since18       = RT13V130,
    Amt_Home_Since18         = RT13V131,
    Amt_Other_Financial_Since18 = RT13V132
  ) %>%
  # keep only child records (type 1), drop parent records (type 2)
  filter(Record_Type == 1, Child_1968_ID < 3000) %>% 
  mutate(
    # Clean education variable
    Child_Education = case_when(
      Child_Education %in% c(98, 99) ~ NA_real_,
      Child_Education == 0           ~ NA_real_,
      TRUE                           ~ as.numeric(Child_Education)
    ),
    
    # Clean age
    Child_Age = case_when(
      Child_Age %in% c(998, 999) ~ NA_real_,
      Child_Age == 0             ~ NA_real_,
      TRUE                       ~ as.numeric(Child_Age)
    ),
    
    # Clean school transfer amount
    # 0 = inapplicable (did not give), 999999998 = DK, 999999999 = NA/refused
    Help_School_Indicator = case_when(
      Amt_School_Since18 > 0 & Amt_School_Since18 < 999999997 ~ 1L,
      Amt_School_Since18 == 0                                  ~ 0L,
      TRUE                                                     ~ NA_integer_
    ),
    Help_School_Amount = case_when(
      Amt_School_Since18 > 0 & Amt_School_Since18 < 999999997 ~ as.numeric(Amt_School_Since18),
      Amt_School_Since18 == 0                                  ~ 0,
      TRUE                                                     ~ NA_real_
    ),
    
    # Clean home transfer amount (for additional analysis)
    Help_Home_Amount = case_when(
      Amt_Home_Since18 > 0 & Amt_Home_Since18 < 999999997 ~ as.numeric(Amt_Home_Since18),
      Amt_Home_Since18 == 0                                ~ 0,
      TRUE                                                 ~ NA_real_
    ),
    
    # Clean other financial help amount
    Help_Other_Amount = case_when(
      Amt_Other_Financial_Since18 > 0 & Amt_Other_Financial_Since18 < 999999997 ~ as.numeric(Amt_Other_Financial_Since18),
      Amt_Other_Financial_Since18 == 0                                          ~ 0,
      TRUE                                                                       ~ NA_real_
    ),
    
    # Annual money transfers (2012) — need to annualize using per-unit codes
    # Per codes: 3=week, 4=two weeks, 5=month, 6=year, 7=other (already annualized by PSID)
    Amt_Parent_Gave_Annual = case_when(
      Amt_Parent_Gave_Money >= 999999997                              ~ NA_real_,
      Amt_Parent_Gave_Money == 0                                      ~ 0,
      Amt_Parent_Gave_Per == 3                                        ~ Amt_Parent_Gave_Money * 52,
      Amt_Parent_Gave_Per == 4                                        ~ Amt_Parent_Gave_Money * 26,
      Amt_Parent_Gave_Per == 5                                        ~ Amt_Parent_Gave_Money * 12,
      Amt_Parent_Gave_Per %in% c(6, 7)                                ~ as.numeric(Amt_Parent_Gave_Money),
      TRUE                                                            ~ NA_real_
    ),
    
    # Degree type from completed education
    # 12 = HS/GED, 13 = <1yr college, 14 = 1yr, 15 = 2yr, 16 = 3yr, 
    # 17 = 4yr college, 18 = 5yr+, 19 = college grad, 20 = post-grad
    degree_type_rt = case_when(
      Child_Education <= 12                          ~ "no_college",
      Child_Education %in% 13:15                     ~ "2yr",     # some college / 2yr
      Child_Education >= 16                          ~ "4yr",     # 4yr degree or more
      TRUE                                           ~ NA_character_
    ),
    
    # Child living arrangement
    child_coresident = as.integer(Child_In_FU == 1)
  ) 

#################################################################################################################################
# Merge and Export
#################################################################################################################################

# Step 1: merge enrollment timing on child
psid_rt13_with_enrollment <- psid_rt13_clean %>%
  left_join(college_enrollment, 
            by = c("Child_1968_ID" = "ER30001", "Child_Person_Number" = "ER30002"))

# Step 2: find closest PSID wave to enrollment start for parent
parent_income_at_decision <- psid_rt13_with_enrollment %>%
  filter(!is.na(enroll_start_year)) %>%
  select(Head_1968_ID, Head_Person_Number, Child_Person_Number, enroll_start_year) %>%
  inner_join(psid_income,
             by = c("Head_1968_ID" = "ER30001", "Head_Person_Number" = "ER30002"),
             relationship = "many-to-many") %>%
  mutate(year_gap = abs(Year - enroll_start_year)) %>%
  group_by(Head_1968_ID, Head_Person_Number, Child_Person_Number) %>%
  slice_min(year_gap, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 3: merge back
psid_edu <- psid_rt13_with_enrollment %>%
  left_join(parent_income_at_decision,
            by = c("Head_1968_ID", "Head_Person_Number", "Child_Person_Number")) %>% 
  mutate(
    enroll_midpoint = round((enroll_start_year + coalesce(enroll_end_year, 2012)) / 2)
  ) %>% 
  left_join(cpi_data, by = ("enroll_midpoint" = "Year")) %>% 
  mutate(
    # expected degree length
    degree_years = case_when(
      degree_type == "2yr" ~ 2,
      degree_type == "4yr" ~ 4,
      TRUE                 ~ NA_real_ ),
    
    # still enrolled at time of RT13
    years_enrolled_so_far = pmax(2012 - enroll_start_year + 1, 1),
    still_enrolled = !enrollment_unobserved &
      (enroll_end_year > 2012 |
         (is.na(enroll_end_year) & !is.na(enroll_start_year) & enroll_start_year <= 2012)),
    
    # tuition help
    scale_factor = case_when(
      still_enrolled & years_enrolled_so_far < degree_years
      ~ pmin(degree_years / years_enrolled_so_far, 4),
      TRUE ~ 1),
    Help_School_Amount_Real = Help_School_Amount * ratio_2010,
    Help_School_Amount_Adj = Help_School_Amount_Real * scale_factor,

    log_educ_exp = log(Help_School_Amount_Adj + 1))

write.csv(psid_edu , "../data/psid_edu.csv")
