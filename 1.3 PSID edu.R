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
# Children
#######################################################################################

psid_children <- psid_ind %>%
  select(ER30001, ER30002, Family_ID, Survey_Year,Relationship_Head, Age, Completed_Education, Employment_Status_Ind) %>% 
  #left_join(psid_fam, multiple = 'all') %>%
  mutate(Year = Survey_Year - 1) %>%
  filter((Relationship_Head %in% 30:43 & Survey_Year == 1994)| (Relationship_Head %in% 30:39 & Survey_Year != 1994), 
         Employment_Status_Ind == 7) %>%
  mutate(
    Age = case_when(
      Age %in% c(99, 999) ~ NA_real_,
      Age == 0             ~ NA_real_,
      TRUE                 ~ as.numeric(Age)
    ),
    # child education completion
    child_college = case_when(
      Completed_Education >= 16 ~ "College Degree",      # 16+ years = BA or more
      Completed_Education >= 13 ~ "Some College",         # 13-15 years = some college
      Completed_Education <= 12 ~ "No College",           # 12 or less = HS or less
      TRUE            ~ NA_character_
    )) %>%
  group_by(Family_ID, Survey_Year) %>%
  summarize(
    n_children_college_age  = sum(Age >= 18 & Age <= 24, na.rm = TRUE),
    n_children_under18      = sum(Age < 18, na.rm = TRUE),
    oldest_child_age        = suppressWarnings(max(Age, na.rm = TRUE)),
    oldest_child_age        = if_else(
      is.infinite(oldest_child_age), NA_real_, oldest_child_age
    ),
    n_children_college      = sum(Completed_Education > 12 & Completed_Education < 17, na.rm = TRUE),
    .groups = "drop"
  )

psid_educ = 
