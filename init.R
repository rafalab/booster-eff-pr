library(tidyverse)
library(lubridate)
library(splines)
library(kableExtra)
library(RColorBrewer)

wait <- 14

min_match_score <- 0.4

vax_status <- brewer.pal(9, "Set1")## colors for vaccination status

plot_manu_levels <- c("UNV", "MOD", "PFR", "JSN", "before", "after")

manu_colors <- c(UNV = vax_status[1], MOD = "#00BFC4", PFR = "#C77CFF", JSN = "#7CAE00", before = hue_pal()(9)[[1]], after = hue_pal()(9)[[6]])

manu_levels <- c("UNV", "MOD", "PFR", "JSN")

manu_labels <- c(UNV = "Unvaccinated", MOD = "mRNA-1273",
                 PFR = "BNT162b2", JSN = "Ad26.COV2.S")

first_day <- make_date(2020, 12, 15)

last_day <- make_date(2021, 10, 15) 

delta_date <- make_date(2021, 6, 15)

all_dates <- tibble(date = seq(first_day, last_day, "days"))

age_starts <- c(12, 18, 25, seq(35, 85, 10))
age_ends <- c(age_starts[-1]-1, Inf)
age_levels <- paste(age_starts, age_ends, sep = "-")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")