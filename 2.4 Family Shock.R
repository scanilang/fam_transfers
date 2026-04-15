library(tidyverse)

psid_clean <- read.csv("../data/psid_clean.csv")

family_shock_table <- psid_clean %>%
  filter(
    Age_recode >= 25, Age_recode <= 45,   # window around 35
    Marital_Status %in% c("Married", "Single"),
    Race_Head %in% c("White", "Black"),
    !is.na(family_type),
    !is.na(Family_Unit_Size)
  ) %>%
  mutate(
    R = if_else(Race_Head == "White", 1L, 2L),
    m = if_else(Marital_Status == "Married", 2L, 1L),
    n = pmin(as.integer(Family_Unit_Size), 5L),
    n = if_else(m == 1L, 1L, n),
    e = case_when(
      Head_College == "No College" ~ 0L,
      TRUE                         ~ 1L   # any college
    ),
    t = family_type,
    age_dist = abs(Age_recode - 35)
  ) %>%
  # One observation per person, closest to age 35
  group_by(ER30001, ER30002) %>%
  slice_min(age_dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(R, e, m, n, t) %>%
  summarize(n_obs = n(), .groups = "drop") %>%
  group_by(R, e) %>%
  mutate(prob = n_obs / sum(n_obs)) %>%
  ungroup()

# Check cell sizes
family_shock_table %>%
  group_by(R, e) %>%
  summarize(
    total_people = sum(n_obs),
    min_cell     = min(n_obs),
    n_cells      = n(),
    .groups = "drop"
  )

write.csv(family_shock_table, "../data/family_shock_table.csv", row.names = FALSE)
