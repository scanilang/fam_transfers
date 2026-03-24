library(tidyverse)

#######################################################################################
# Read in raw data 
#######################################################################################
tas_raw <- read.csv("./data/transition_adulthood.csv")
cpi_data <- openxlsx::read.xlsx("./data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 

#######################################################################################
# Clean data 
#######################################################################################
tas_clean <- tas_raw %>% 
  mutate(Year = Survey_Year - 1) %>% 
  left_join(cpi_data %>% select(year, ratio_2010), by = c("Year" = "year")) %>% 
  mutate(across(matches("^Help_.*Amount_Parents$"), ~ na_if(na_if(., 9999998), 9999999)),
         across(matches("^Help_.*Amount_Parents$"), ~ . * ratio_2010),
         Tuition_Amount = Tuition_Amount * ratio_2010,
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Gifts_Inheritance"), ~ na_if(na_if(., 9999998), 9999999)),
         across(starts_with("Tuition_Amount"), ~ na_if(na_if(., 9999998), 9999999))) 


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
# Help with tuition
#######################################################################################
tas_clean_college = tas_clean %>% 
  filter(Race %in% c(1,2),  Race_2nd == 0) %>% 
  filter(!is.na(Help_Tuition_Amount_Parents), !is.na(Tuition_Amount)) %>% 
  mutate(college_status = case_when(Enrollment_Status == 9 ~ "Enrolled",
                                    Enrollment_Status == 11 ~ "Enrolled, college degree",
                                    Enrollment_Status %in% 6:7 ~ "College degree",
                                    Enrollment_Status %in% 4:5 ~ "Some college",
                                    Enrollment_Status %in% 1:3 ~ "No college")) %>% 
  group_by(Race, college_status) %>% 
  summarize(n = n(),
            avg_help = mean(Help_Tuition_Amount_Parents),
            avg_tution = mean(Tuition_Amount))
