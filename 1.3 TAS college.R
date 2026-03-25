library(tidyverse)

#######################################################################################
# Read in raw data 
#######################################################################################
tas_raw <- read.csv("../data/transition_adulthood.csv")
cpi_data <- openxlsx::read.xlsx("../data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 
psid_clean <- read.csv("../data/psid_clean.csv")
psid_hh <- read.csv("../data/famtransfer_famfile.csv")

#######################################################################################
# Clean data 
#######################################################################################
tas_clean <- tas_raw %>% 
  mutate(Year = Survey_Year - 1) %>% 
  left_join(cpi_data %>% select(year, ratio_2010), by = c("Year" = "year")) %>% 
  mutate(across(matches("^Help_.*Amount_Parents$"), ~ na_if(na_if(., 9999998), 9999999)),
         across(matches("^Help_.*Percent_Parents$"), ~ na_if(na_if(., 998), 999)),
         across(matches("^Help_.*Amount_Parents$"), ~ . * ratio_2010),
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Tuition_Amount"), ~ na_if(na_if(., 9999998), 9999999)),
         Tuition_Amount = Tuition_Amount * ratio_2010,
         Student_Loan_Amount = Student_Loan_Amount * ratio_2010) %>%
  mutate(full_year_4yr = case_when(Degree_Type_Seeking_1 == 2 & (First_Year_Attended_1 < Year | First_Year_Attended_1 == Year & First_Month_Attended_1 == 1) & 
                                     (Last_Year_Attended_1 > Year | Last_Year_Attended_1 == 0 | Last_Year_Attended_1 == Year & Last_Month_Attended_1 == 12) ~ 1,
                                   Degree_Type_Seeking_2 == 2 & (First_Year_Attended_2 < Year | First_Year_Attended_2 == Year & First_Month_Attended_2 == 1) & 
                                     (Last_Year_Attended_2 > Year | Last_Year_Attended_2 == 0 | Last_Year_Attended_2 == Year & Last_Month_Attended_2 == 12) ~ 1,
                                   Degree_Type_Seeking_1 == 2 & Degree_Type_Seeking_2 == 2 & First_Year_Attended_1 < Year & Last_Year_Attended_1 ==Year & 
                                     (Last_Year_Attended_2 == 0 | Last_Year_Attended_2 > Year | Last_Year_Attended_2 == Year & Last_Month_Attended_2 == 12) ~ 1,
                                   TRUE ~ 0),
         full_year_2yr = case_when(Degree_Type_Seeking_1 == 1 & (First_Year_Attended_1 < Year | First_Year_Attended_1 == Year & First_Month_Attended_1 == 1) & 
                                     (Last_Year_Attended_1 > Year | Last_Year_Attended_1 == 0 | Last_Year_Attended_1 == Year & Last_Month_Attended_1 == 12) ~ 1,
                                   Degree_Type_Seeking_2 == 1 & (First_Year_Attended_2 < Year | First_Year_Attended_2 == Year & First_Month_Attended_2 == 1) & 
                                     (Last_Year_Attended_2 > Year | Last_Year_Attended_2 == 0 | Last_Year_Attended_2 == Year & Last_Month_Attended_2 == 12) ~ 1,
                                   Degree_Type_Seeking_1 == 1 & Degree_Type_Seeking_2 == 1 & First_Year_Attended_1 < Year & Last_Year_Attended_1 ==Year & 
                                     (Last_Year_Attended_2 == 0 | Last_Year_Attended_2 > Year | Last_Year_Attended_2 == Year & Last_Month_Attended_2 == 12) ~ 1,
                                   TRUE ~ 0)) %>% 
  left_join(psid_clean %>% select(Family_ID, Survey_Year, Total_Income_Head_Spouse, Total_Family_Income, Family_Unit_Size))


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

#######################################################################################
# Tuition Amount
#######################################################################################
tuition_amount = tas_clean %>% 
  select(Family_ID, ID_number, Survey_Year, Year, Enrollment_Status, Part_or_Full_Time_Student, Degree_Type_Seeking_1, Degree_Type_Seeking_2, 
         First_Month_Attended_1, First_Year_Attended_1, Last_Month_Attended_1, Last_Year_Attended_1,
         First_Month_Attended_2, First_Year_Attended_2, Last_Month_Attended_2, Last_Year_Attended_2, Tuition_Amount, full_year_4yr, full_year_2yr) %>% 
  # Full Time Student and Full Year Enrolled
  filter(Part_or_Full_Time_Student == 1) %>% 
  group_by(Degree_Type_Seeking_1) %>% 
  dplyr::summarize(avg_tution)
            
#######################################################################################
# Help with tuition
#######################################################################################
tas_clean_college = tas_clean %>% 
  arrange(Family_ID, ID_number) %>% 
  select(Family_ID, ID_number, Survey_Year, Race, Race_2nd, Enrollment_Status,Tuition_Amount, Help_Tuition, 
         Help_Tuition_Amount_Parents, Help_Tuition_Percent_Parents, Student_Loan_Amount, Help_Student_Loan, 
         Help_Student_Loan_Amount_Parents, Help_Student_Loan_Percent_Parents,
         Total_Income_Head_Spouse, Total_Family_Income, Family_Unit_Size) %>% 
  filter(Race %in% c(1,2),  Race_2nd == 0) %>% 
  filter(!is.na(Help_Tuition_Amount_Parents), !is.na(Tuition_Amount)) %>% 
  mutate(college_status = case_when(Enrollment_Status == 9 ~ "Enrolled, no prior degree",
                                    Enrollment_Status == 10 ~ "Enrolled, prior degree",
                                    Enrollment_Status == 11 ~ "Enrolled, college degree",
                                    Enrollment_Status %in% 6:7 ~ "College degree",
                                    Enrollment_Status %in% 4:5 ~ "Some college",
                                    Enrollment_Status %in% c(1:3,95) ~ "No college"),
         Help_Tuition_Amount_Parents = if_else(Help_Tuition_Amount_Parents != 0, Help_Tuition_Amount_Parents,
                                               Tuition_Amount * (Help_Tuition_Percent_Parents/100))) 

avgs = tas_clean_college %>% 
  group_by(Race, college_status) %>% 
  summarize(n = n(),
            avg_help = mean(Help_Tuition_Amount_Parents),
            avg_tution = mean(Tuition_Amount))


# Probability of receiving help with tuition
tuition_prob = glm(Help_Tuition ~ Race + Total_Income_Head_Spouse + Family_Unit_Size, data = tas_clean_college)

# Amount of help with tuition
tuition_help = lm(Help_Tuition_Amount_Parents ~ Race + Total_Income_Head_Spouse + Family_Unit_Size, data = tas_clean_college)
