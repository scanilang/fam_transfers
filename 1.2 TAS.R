library(tidyverse)

#######################################################################################
# Read in raw data 
#######################################################################################
tas_raw <- read.csv("../data/transition_adulthood.csv")
cpi_data <- openxlsx::read.xlsx("../data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 
psid_clean <- read.csv("../data/psid_clean.csv") 

#######################################################################################
# Clean data 
#######################################################################################
tas_clean <- tas_raw %>% 
  mutate(Year = Survey_Year - 1) %>% 
  left_join(cpi_data %>% select(year, ratio_2010), by = c("Year" = "year")) %>% 
  mutate(across(matches("^Help_.*Amount_Parents$"), ~ na_if(na_if(., 9999998), 9999999)),
         across(ends_with("_Amount"), ~ na_if(na_if(., 9999998), 9999999)),
         across(matches("Help_Personal_Loan_Amount"), ~ na_if(na_if(., 9999998), 9999999)),
         across(matches("Help_Other_Amount"), ~ na_if(na_if(., 9999998), 9999999)),
         across(matches("^Help_.*Percent_Parents$"), ~ na_if(na_if(., 998), 999)),
         across(matches("^Help_.*Amount_Parents$"), ~ . * ratio_2010),
         across(ends_with("_Amount"), ~ . * ratio_2010),
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Tuition_Amount"), ~ na_if(na_if(., 9999998), 9999999)),
         Help_Tuition = if_else(Help_Tuition == 1, 1, 0)) %>%
  mutate(college_status = case_when(Enrollment_Status == 9 ~ "Enrolled, no prior degree",
                                    Enrollment_Status == 10 ~ "Enrolled, prior degree",
                                    Enrollment_Status == 11 ~ "Enrolled, college degree",
                                    Enrollment_Status %in% 6:7 ~ "College degree",
                                    Enrollment_Status %in% 4:5 ~ "Some college",
                                    Enrollment_Status %in% c(1:3,95) ~ "No college"),
         Help_Tuition_Amount_Parents = if_else(Help_Tuition_Amount_Parents != 0, Help_Tuition_Amount_Parents,
                                               Tuition_Amount * (Help_Tuition_Percent_Parents/100))) %>% 
  mutate(
    # active during reference year
    active_in_year_1 = (
      First_Year_Attended_1 <= Year &
        (Last_Year_Attended_1 >= Year | Last_Year_Attended_1 == 0)
    ),
    active_in_year_2 = (
      First_Year_Attended_2 <= Year &
        (Last_Year_Attended_2 >= Year | Last_Year_Attended_2 == 0)
    ),
    
    # degree type — loose condition, just needs to be active during year
    degree_type = case_when(
      active_in_year_1 & Degree_Type_Seeking_1 == 2 ~ "4yr",
      active_in_year_1 & Degree_Type_Seeking_1 == 1 ~ "2yr",
      active_in_year_2 & Degree_Type_Seeking_2 == 2 ~ "4yr",
      active_in_year_2 & Degree_Type_Seeking_2 == 1 ~ "2yr",
      TRUE ~ NA_character_
    ),
    
    transfer_contiguous = (
      active_in_year_1 & active_in_year_2 &
        (
          (First_Year_Attended_2 == Last_Year_Attended_1 &
             First_Month_Attended_2 <= Last_Month_Attended_1 + 2) |
            (First_Year_Attended_2 == Last_Year_Attended_1 + 1 &
               Last_Month_Attended_1 == 12 &
               First_Month_Attended_2 == 1)
        )
    ),
    
    # full year flags — strict condition
    enrolled_full_start_1 = (
      First_Year_Attended_1 < Year |
        (First_Year_Attended_1 == Year & First_Month_Attended_1 == 1)
    ),
    enrolled_full_start_2 = (
      First_Year_Attended_2 < Year |
        (First_Year_Attended_2 == Year & First_Month_Attended_2 == 1)
    ),
    enrolled_full_end_1 = (
      Last_Year_Attended_1 > Year |
        Last_Year_Attended_1 == 0   |
        (Last_Year_Attended_1 == Year & Last_Month_Attended_1 == 12)
    ),
    enrolled_full_end_2 = (
      Last_Year_Attended_2 > Year |
        Last_Year_Attended_2 == 0   |
        (Last_Year_Attended_2 == Year & Last_Month_Attended_2 == 12)
    ),
    
    full_year = (
      (active_in_year_1 & enrolled_full_start_1 & enrolled_full_end_1) |
        (active_in_year_2 & enrolled_full_start_2 & enrolled_full_end_2) |
        # transfer case: same degree type, contiguous
        (active_in_year_1 & active_in_year_2 &
           Degree_Type_Seeking_1 == Degree_Type_Seeking_2 &
           enrolled_full_start_1 & enrolled_full_end_2 &
           transfer_contiguous)
    ),
    
    # one semester: active during year but not full year
    # either started mid-year or graduated mid-year
    one_semester = (
      !is.na(degree_type) &
        !full_year &
        (
          # started fall of reference year, still going
          (active_in_year_1 & First_Year_Attended_1 == Year &
             First_Month_Attended_1 %in% c(8, 9) &
             (Last_Year_Attended_1 == 0 | Last_Year_Attended_1 > Year)) |
            (active_in_year_2 & First_Year_Attended_2 == Year &
               First_Month_Attended_2 %in% c(8, 9) &
               (Last_Year_Attended_2 == 0 | Last_Year_Attended_2 > Year)) |
            # graduated spring of reference year
            (active_in_year_1 & enrolled_full_start_1 &
               Last_Year_Attended_1 == Year &
               Last_Month_Attended_1 %in% c(5, 6)) |
            (active_in_year_2 & enrolled_full_start_2 &
               Last_Year_Attended_2 == Year &
               Last_Month_Attended_2 %in% c(5, 6))
        )
    ),
    
    # annualized tuition
    tuition_annual = case_when(
      !is.na(degree_type) & full_year    ~ Tuition_Amount,
      !is.na(degree_type) & one_semester ~ Tuition_Amount * 2,
      TRUE ~ NA_real_
    )) %>% 
  mutate(
      # financial constraint — cant afford
      cant_afford = as.integer(
        Reason_Not_In_College_1 == 8 |
          Reason_Not_In_College_2 == 8 |
          Reason_Not_In_College_3 == 8),
      # opportunity cost constraint
      needed_to_work = as.integer(
        Reason_Not_In_College_1 == 9 |
          Reason_Not_In_College_2 == 9 |
          Reason_Not_In_College_3 == 9),
      # any financial/economic constraint
      financially_constrained = as.integer(
        cant_afford == 1 | needed_to_work == 1)) %>% 
  mutate(Year_bins_edu = case_when(
    Year %in% 2004:2010 ~ "2004-2010",
    Year %in% 2011:2020 ~ "2011-2020",
    TRUE ~ NA_character_
  )) %>% 
  left_join(psid_clean %>% rename(Marital_Status_Parents = Marital_Status) %>% 
              select(Family_ID, Survey_Year, Race_Head, Head_College, log_asset_income, log_nonasset_income, Marital_Status_Parents, Family_Unit_Size, family_type),
              by = c("Family_ID" = "Family_ID", "Survey_Year" = "Survey_Year"))

#######################################################################################
# Tuition Amount
#######################################################################################
tuition_amount = tas_clean %>% 
  select(Family_ID, ID_number, Survey_Year, Individual_Weight, Year, Enrollment_Status, college_status, Part_or_Full_Time_Student, Degree_Type_Seeking_1, Degree_Type_Seeking_2, 
         First_Month_Attended_1, First_Year_Attended_1, Last_Month_Attended_1, Last_Year_Attended_1,
         First_Month_Attended_2, First_Year_Attended_2, Last_Month_Attended_2, Last_Year_Attended_2, Tuition_Amount,
         degree_type, tuition_annual) %>% 
  # Full Time Student and Full Year Enrolled
  filter(Part_or_Full_Time_Student == 1,
         !is.na(Individual_Weight),
         tuition_annual > 0,
         college_status != "No college") %>% 
  group_by(degree_type) %>% 
  summarize(
    avg_tuition = weighted.mean(tuition_annual, w = Individual_Weight, na.rm = TRUE),
    n = n(),
    .groups = "drop")

#######################################################################################
# Living Support During School (non-tuition transfers)
#######################################################################################

tas_school_transfers <- tas_clean %>%
  filter(Enrollment_Status %in% c(9, 10),      # currently enrolled
         Part_or_Full_Time_Student == 1, 
         Survey_Year >= 2017,
         Race_Head %in% c("White", "Black")) %>%     # full time
  mutate(
    # Construct parental amount for each category
    # Use amount if available, otherwise pct * total
    Help_Rent_Mortgage_Parents = case_when(
      Help_Rent_Mortgage_Amount_Parents > 0 ~ Help_Rent_Mortgage_Amount_Parents,
      Help_Rent_Mortgage_Percent_Parents > 0 ~ Rent_Mortgage_Amount * (Help_Rent_Mortgage_Percent_Parents / 100),
      Help_Rent_Mortgage == 1 ~ NA_real_,  # helped but amount unknown
      TRUE ~ 0
    ),
    
    Help_House_Parents = case_when(
      Help_House_Amount_Parents > 0 ~ Help_House_Amount_Parents,
      Help_House_Percent_Parents > 0 ~ House_Amount * (Help_House_Percent_Parents / 100),
      Help_House == 1 ~ NA_real_,  # helped but amount unknown
      TRUE ~ 0
    ),
    
    Help_Bills_Parents = case_when(
      Help_Bills_Amount_Parents > 0 ~ Help_Bills_Amount_Parents,
      Help_Bills_Percent_Parents > 0 ~ Bills_Amount * (Help_Bills_Percent_Parents / 100),
      Help_Bills == 1 ~ NA_real_,
      TRUE ~ 0
    ),
    Help_Vehicle_Parents = case_when(
      Help_Personal_Vehicle_Amount_Parents > 0 ~ Help_Personal_Vehicle_Amount_Parents,
      Help_Personal_Vehicle_Percent_Parents > 0 ~ Personal_Vehicle_Amount * (Help_Personal_Vehicle_Percent_Parents / 100),
      Help_Personal_Vehicle == 1 ~ NA_real_,
      TRUE ~ 0
    ),
    Help_Other_Parents = case_when(
      Help_Other_Amount > 0 & Help_Other == 1 ~ Help_Other_Amount,
      Help_Other == 1 ~ NA_real_,
      TRUE ~ 0
    ),
    Help_Personal_Loan_Parents =  case_when(
      Help_Personal_Loan_Amount > 0 ~ Help_Personal_Loan_Amount,
      Help_Personal_Loan == 1 ~ NA_real_,
      TRUE ~ 0
    ),
    
    # Total non-tuition living support
    Living_Support = Help_Rent_Mortgage_Parents + Help_Bills_Parents + Help_Vehicle_Parents + Help_Other_Parents + Help_Personal_Loan_Parents + Help_House_Parents,
  
    Living_Support_Indicator = as.integer(Living_Support > 0),
    
    # Also compute total support including tuition for comparison
   # Total_Support = Living_Support 
    
  ) %>% 
  select(Family_ID, Survey_Year, Race, Race_Head, Head_College, log_asset_income, log_nonasset_income, Marital_Status_Parents, Family_Unit_Size, family_type, 
         Living_Support_Indicator, Help_Rent_Mortgage_Parents,Help_Bills_Amount_Parents, Help_Vehicle_Parents, Help_Personal_Loan_Parents, 
         Help_Other_Parents,Help_House_Parents, Individual_Weight, Living_Support )

write.csv(tas_school_transfers, "../data/tas_school_transfers.csv")

# Descriptives
tas_school_transfers %>%
  filter(Race_Head %in% c(1, 2)) %>%
  mutate(Race_label = if_else(Race_Head == 1, "White", "Black")) %>%
  group_by(Race_label) %>%
  summarize(
    n = n(),
    share_any_living = weighted.mean(Living_Support_Indicator, Individual_Weight, na.rm = TRUE),
    mean_living_if_pos = weighted.mean(Living_Support[Living_Support > 0], 
                                        Individual_Weight[Living_Support > 0], na.rm = TRUE),
     mean_total_support = weighted.mean(Total_Support, Individual_Weight, na.rm = TRUE),
    .groups = "drop"
  )
  
#######################################################################################
# Some Summaries
#######################################################################################



tas_clean_summary = tas_clean %>% 
  filter(Race %in% c(1,2)) %>% 
  group_by(Race) %>% 
  summarize(n = n(),
            received_personal_laon = sum(Help_Personal_Loan == 1) / n,
            receoved_help_rent_mortgage = sum(Help_Rent_Mortgage == 1) / n,
            received_help_bills = sum(Help_Bills == 1) /n )


tas_clean_college = tas_clean %>% 
  filter(Race %in% c(1,2)) %>% 
  filter(Enrollment_Status >= 3) %>% 
  filter(Enrollment_Status < 99) %>%
  group_by(Race) %>% 
  summarize(n = n(),
            received_help_tuition = sum(Help_Tuition == 1) / n)

tas_clean_college = tas_clean %>% 
  filter(Race %in% c(1,2), Race_2nd == 0) %>% 
  filter(Help_Tuition_Amount_Parents > 0) %>% 
  group_by(Race) %>% 
  summarize(n = n(),
            received_help_tuition = mean(Help_Tuition_Amount_Parents))


# empirical moment from TAS
cant_afford_moment <- tas_clean %>%
  filter(Race %in% c(1, 2),
         Enrollment_Status %in% c(1, 2, 3, 4, 5, 95)) %>%  # not enrolled
  group_by(Race) %>%
  summarize(
    pct_cant_afford = weighted.mean(cant_afford, Individual_Weight, na.rm = TRUE)
  )


