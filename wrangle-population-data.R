library(tidyverse)
library(lubridate)
library(data.table)

load("rdas/dates.rda")
pr_pop <- 3285874 ## population of puerto rico

age_starts <- as.numeric(str_extract(age_levels, "\\d+"))
breaks <- sort(age_starts)
labels <- c(paste(breaks[-length(breaks)], c(breaks[-1]-1), sep="-"),
            paste0(breaks[length(breaks)], "+"))
pop_by_age_gender  <-
    readxl::read_xlsx("data/prc-est2019-syasex.xlsx",  
                               skip = 4) %>%
  select("...1", "Male...36", "Female...37")  %>% 
  setNames(c("ageRange", "M", "F")) %>%
  mutate(ageRange= str_remove(ageRange, "."))  %>% 
  filter(str_detect(ageRange, "^\\d")) %>% 
  mutate(ageRange = str_remove(ageRange,"\\+")) %>%
  mutate(ageRange = as.numeric(ageRange)) %>% 
  pivot_longer(-ageRange, names_to="gender", values_to = "poblacion") 

pr_adult_pop <- pop_by_age_gender %>% filter(ageRange>= 18) %>%
  summarize(total = sum(poblacion)) %>%
  pull(total)
old_pop <- sum(pop_by_age_gender$poblacion)

pop_by_age_gender <- pop_by_age_gender %>%
  mutate(poblacion = poblacion*pr_pop/old_pop) %>% 
  mutate(ageRange = factor(as.character(cut(ageRange, c(breaks, Inf), right = FALSE,
                                            labels = labels)), levels = labels)) %>%
  filter(!is.na(ageRange)) %>%
  group_by(ageRange, gender) %>%
  summarize(poblacion = sum(poblacion), .groups = "drop") %>%
  as.data.table() 

save(pr_pop, pr_adult_pop, pop_by_age_gender, file = "rdas/population-tabs.rda")


