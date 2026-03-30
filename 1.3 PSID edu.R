library(dplyr)

#######################################################################################
# Read in raw data 
#######################################################################################
psid_fam<- read.csv("../data/famtransfer_famfile.csv") %>% rename(ER30001 = ID_1968)
psid_ind<- read.csv("../data/famtransfer_indfile.csv") %>% filter(ER30001 < 3000) %>% select(-Interview_Number)
cpi_data <- openxlsx::read.xlsx("../data/bls_CPI.xlsx") %>% 
  mutate(ratio_2010 = 218.056 / annual.avg) 

# cleaned psid data
psid_clean <- read.csv("../data/psid_clean.csv")

# function to take the mode
getmode <- function(v) {
  v = na.omit(v)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#######################################################################################
# Children in College
#######################################################################################
psid_children <- psid_ind %>%
  select(ER30001, ER30002, Family_ID, Survey_Year,
         Relationship_Head, Age, Completed_Education, 
         Employment_Status_Ind) %>% 
  mutate(Year = Survey_Year - 1) %>%
  mutate(
    Age = case_when(
      Age %in% c(99, 999) ~ NA_real_,
      Age == 0             ~ NA_real_,
      TRUE                 ~ as.numeric(Age)
    ),
    Completed_Education = if_else(
      Completed_Education %in% c(98, 99), NA_real_,
      as.numeric(Completed_Education)
    )
  ) %>%
  # get eventual education for each child
  arrange(ER30001, ER30002, Survey_Year) %>%
  group_by(ER30001, ER30002) %>%
  mutate(
    max_education = max(Completed_Education, na.rm = TRUE),
    max_education = if_else(
      is.infinite(max_education), NA_real_, max_education),
    
    # first age where completed education reaches 2yr threshold (14 years)
    reached_2yr       = Completed_Education >= 14,
    age_at_2yr        = if_else(reached_2yr, Age, NA_real_),
    age_graduated_2yr = min(age_at_2yr, na.rm = TRUE),
    age_graduated_2yr = if_else(
      is.infinite(age_graduated_2yr), NA_real_, age_graduated_2yr),
    
    # first age where completed education reaches 4yr threshold (16 years)
    reached_4yr       = Completed_Education >= 16,
    age_at_4yr        = if_else(reached_4yr, Age, NA_real_),
    age_graduated_4yr = min(age_at_4yr, na.rm = TRUE),
    age_graduated_4yr = if_else(
      is.infinite(age_graduated_4yr), NA_real_, age_graduated_4yr),
    
    # graduation age based on eventual degree type
    age_graduated = case_when(
      max_education %in% 13:14 ~ age_graduated_2yr,
      max_education >= 15      ~ age_graduated_4yr,
      TRUE                     ~ NA_real_
    ),
    
    # enrollment age — backtrack from graduation
    age_enrolled_est = case_when(
      max_education %in% 13:14 ~ age_graduated_2yr - 2,
      max_education >= 15      ~ age_graduated_4yr - 4,
      TRUE                     ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  # filter to children in household
  filter(
    (Relationship_Head %in% 30:43 & Survey_Year == 1994) | 
      (Relationship_Head %in% 30:39 & Survey_Year != 1994)
  ) %>%
  # summarize household level
  group_by(Family_ID, Survey_Year) %>%
  summarize(
    n_children_total       = n(),
    n_children_under18     = sum(Age < 18, na.rm = TRUE),
    n_children_college_age = sum(Age >= 18 & Age <= 26, na.rm = TRUE),
    oldest_child_age       = suppressWarnings(max(Age, na.rm = TRUE)),
    oldest_child_age       = if_else(
      is.infinite(oldest_child_age), NA_real_, oldest_child_age
    ),
    n_children_college     = sum(
      Age >= 18 & Age <= 26 &
        Completed_Education > 12 & Completed_Education <= 16,
      na.rm = TRUE
    ),
    any_child_in_college   = as.integer(n_children_college > 0),
    # use max education to identify degree type
    # child eventually completing 13-14 years → 2yr degree
    # child eventually completing 15-16 years → 4yr degree
    n_children_2yr         = sum(
      Age >= 18 & Age <= 26 &
        Completed_Education > 12 &
        max_education %in% 13:14,
      na.rm = TRUE
    ),
    n_children_4yr         = sum(
      Age >= 18 & Age <= 26 &
        Completed_Education > 12 &
        max_education >= 15,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

#######################################################################################
# Expenditure per Child in College
#######################################################################################

psid_educ <- psid_clean %>%
  left_join(psid_children, by = c("Family_ID", "Survey_Year")) %>%
  filter(Race_Head %in% c("White", "Black")) %>%
  arrange(ER30001, ER30002, Survey_Year) %>%
  group_by(ER30001, ER30002) %>%
  mutate(
    # flag when first child enters college
    child_in_college      = if_else(
      is.na(any_child_in_college), 0L, any_child_in_college),
    child_in_college_lag  = lag(child_in_college),
    
    # transition into college: 0 -> 1
    college_entry         = as.integer(
      child_in_college == 1 & 
        (is.na(child_in_college_lag) | child_in_college_lag == 0)
    ),
    
    # change in education expenditure
    educ_exp_lag          = lag(Total_Education_Expenditure),
    delta_educ_exp        = Total_Education_Expenditure - 
      replace_na(educ_exp_lag, 0),
    
    # number of children in college this period vs last
    n_college_lag         = lag(n_children_college),
    delta_n_college       = n_children_college - 
      replace_na(n_college_lag, 0)
  ) %>%
  ungroup()

# check jump
psid_educ %>%
  filter(Race_Head %in% c("White", "Black")) %>%
  mutate(
    college_status = case_when(
      child_in_college == 0 & child_in_college_lag == 0 ~ "no college",
      child_in_college == 1 & child_in_college_lag == 0 ~ "entry year",
      child_in_college == 1 & child_in_college_lag == 1 ~ "enrolled",
      child_in_college == 0 & child_in_college_lag == 1 ~ "exit year",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Race_Head, college_status) %>%
  summarize(
    n           = n(),
    avg_educ    = weighted.mean(
      Total_Education_Expenditure, na.rm = TRUE),
    avg_delta   = weighted.mean(
      delta_educ_exp, na.rm = TRUE),
    .groups = "drop"
  )


write.csv(psid_edu , "../data/psid_edu.csv")