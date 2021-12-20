library(data.table)
library(tidyverse)
library(lubridate)
library(scales)

manu_levels <- c("UNV", "MOD", "PFR", "JSN")

load("rdas/population-tabs.rda")
load("rdas/dates.rda")
load("rdas/cases_3_months.rda")
load("rdas/dat_vax.rda")

first_booster_day <- make_date(2021, 8, 13)
first_jnj_booster_day <- make_date(2021, 10, 22)
analysis_last_day <- last_day - weeks(2)
  
##susceptible population
pop_by_age_gender <- cases %>% 
  filter(gender %in% c("F", "M") & !is.na(ageRange)) %>%
  left_join(pop_by_age_gender, by = c("ageRange", "gender")) %>%
  mutate(poblacion = poblacion - cases) %>%
  select(-cases) %>%
  filter(date >= first_day & date <= analysis_last_day)

rm(cases)


dat_vax <- dat_vax[!manu_1 %in% c("ATZ","OTR") & 
          !manu_2 %in% c("ATZ","OTR") &
          !manu_3 %in% c("ATZ","OTR"), ] 
dat_vax$manu_1 <- droplevels(dat_vax$manu_1)
dat_vax$manu_2 <- droplevels(dat_vax$manu_2)
dat_vax$manu_3 <- droplevels(dat_vax$manu_3)

dat_vax[estado %in% c("PR", "No reportado") , c("vax_date", "booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(date_2, date_3, manu_3, insert_date_3, proveedor_3, ageRange_3)]
dat_vax[manu_1 == "JSN", 
        c("vax_date", "booster_date", "ageRange_2", "manu_2", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(date_1, date_2, ageRange_1, manu_1, manu_2, date_2, proveedor_2, ageRange_2)]
dat_vax[!is.na(vax_date), vax_date := vax_date + days(14)]
dat_vax[vax_date > today(), vax_date := NA] ## not fully vax yet

## to early to be a booster
dat_vax[!is.na(booster_date) & manu_1 != "JSN" & (booster_date < first_booster_day | booster_date - date_1 < days(180)), 
        c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]
dat_vax[!is.na(booster_date) & manu_1 == "JSN" & (booster_date < first_jnj_booster_day | booster_date - date_1 < days(60)), 
        c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]
# not vaxed, can't be boosted
dat_vax[!is.na(booster_date) & is.na(vax_date), 
        c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]

## for convinience when using data.table group_by we define date
dat_vax[, date := vax_date]
dat_vax <- dat_vax[(is.na(proveedor_1) | proveedor_1!="Correccional") & 
                     (is.na(proveedor_2) | proveedor_2!="Correccional"),]

### Compute number with and without booster for each day after vaxinated
the_ndays <- as.numeric(analysis_last_day - min(dat_vax$vax_date, na.rm = TRUE))
the_booster_ndays <- as.numeric(analysis_last_day - first_booster_day)

compute_date_comb_counts <- function(tab, ndays = the_ndays){
  the_date <- unique(tab$date)
  all_dates <- data.table(booster_date = seq(the_date, min(the_date + the_ndays, analysis_last_day), "days"))
  ret <- tab[!is.na(booster_date), .(poblacion = .N), keyby = .(booster_date)]
  ret <- merge(all_dates, ret, by = "booster_date", all.x = TRUE)
  ret[is.na(ret)] <- 0
  ret <- ret[, poblacion := cumsum(poblacion)]
  ret[, poblacion := nrow(tab) - poblacion]
  setnames(ret, "booster_date", "date")
  return(ret)
}

pop_vax <- dat_vax[!is.na(vax_date) & vax_date <= analysis_last_day &
                     gender %in% c("F", "M") & !ageRange_1 %in% c("0-11", "12-17"),
                   compute_date_comb_counts(tab = .SD), keyby = c("vax_date", "manu_1", "ageRange_1", "gender")]
setnames(pop_vax, c("vax_date", "manu", "ageRange", "gender", "date", "poblacion"))
setcolorder(pop_vax, c("date", "vax_date", "manu", "ageRange", "gender", "poblacion"))
pop_vax$manu <- factor(pop_vax$manu, levels = manu_levels)
pop_vax$ageRange <- factor(pop_vax$ageRange, levels = age_levels)
pop_vax <- pop_vax[order(manu, ageRange, gender, date, vax_date)]


## Compute unvaxinated
all_dates <- data.table(date = seq(first_day, analysis_last_day, "days"))
all_combs <- CJ(date = all_dates$date, 
                ageRange = levels(dat_vax$ageRange_1)[-c(1,2)],
                gender = levels(dat_vax$gender),
                manu = levels(dat_vax$manu_1))
## First one dose
daily_counts_onedose_age_gender_manu <- dat_vax[!is.na(date_1), .(onedose = .N),
                                                keyby = .(date_1, ageRange_1, gender, manu_1)]
names(daily_counts_onedose_age_gender_manu) <- c("date", "ageRange", "gender", "manu", "onedose")
daily_counts_onedose_age_gender_manu <- merge(all_combs, 
                                              daily_counts_onedose_age_gender_manu, 
                                              all.x=TRUE)
daily_counts_onedose_age_gender_manu[is.na(daily_counts_onedose_age_gender_manu)] <- 0 
counts_onedose_age_gender_manu <- copy(daily_counts_onedose_age_gender_manu)
counts_onedose_age_gender_manu[, onedose := cumsum(onedose), keyby = .(ageRange, gender, manu)]
setnames(counts_onedose_age_gender_manu, "onedose", "n")

## to compute unvax, take population and remove those with at least one dose and those infected last 3 monhts
pop_unvax <- counts_onedose_age_gender_manu[gender %in% c("F", "M") & 
                                              !ageRange %in% c("0-11", "12-17")  &
                                              date <= analysis_last_day,
                                            .(total = sum(n)), 
                                            keyby = .(date, ageRange, gender)]
pop_unvax <- merge(pop_unvax, pop_by_age_gender, all.x = TRUE, 
                   by = c("date", "ageRange", "gender")) 
pop_unvax[, poblacion := poblacion - total]
pop_unvax <- pop_unvax[,c("date", "ageRange","gender", "poblacion")]

### Booster counts
daily_counts_booster_age_gender_manu <- dat_vax[!is.na(booster_date), .(poblacion = .N),
                                                keyby = .(booster_date, ageRange_1, gender, manu_1, booster_manu)]
names(daily_counts_booster_age_gender_manu) <- c("booster_date", "ageRange", "gender", "manu_1", "booster_manu", "poblacion")

## COMPUTE DATA for WAINING

### CASES

load("rdas/dat_cases_vax.rda")
dat_cases <- dat_cases_vax[!manu_1 %in% c("ATZ","OTR") & 
                             !manu_2 %in% c("ATZ","OTR") &
                             !manu_3 %in% c("ATZ","OTR") &
                             (is.na(estado) | estado %in% c("PR", "No reportado")),]
dat_cases <- dat_cases[(is.na(proveedor_1) | proveedor_1 != "Correccional") & 
                         (is.na(proveedor_2) | proveedor_2 != "Correccional"),]
dat_cases$manu_1 <- factor(dat_cases$manu_1, levels = manu_levels)
dat_cases$manu_2 <- factor(dat_cases$manu_2, levels = manu_levels)
dat_cases$manu_3 <- factor(dat_cases$manu_3, levels = manu_levels)


dat_cases[manu_1 != "JSN", c("vax_date", "booster_date", "booster_manu") := .(date_2, date_3, manu_3)]
dat_cases[manu_1 == "JSN", 
          c("vax_date", "booster_date", "manu_2", "booster_manu") := .(date_1, date_2, manu_1, manu_2)] 
dat_cases$manu <- factor(replace_na(as.character(dat_cases$manu_1), "UNV"), levels = manu_levels)
dat_cases[!is.na(vax_date), vax_date := vax_date + days(14)]
dat_cases[!is.na(booster_date) & manu != "JSN" & (booster_date < first_booster_day | booster_date - date_1 < days(180)), 
          c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]
dat_cases[!is.na(booster_date) & manu == "JSN" & (booster_date < first_jnj_booster_day | booster_date - date_1 < days(60)), 
          c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]
# not vaxed, can't be boosted
dat_cases[!is.na(booster_date) & is.na(vax_date), 
        c("booster_date", "booster_manu", "booster_insert_date", "booster_proveedor", "booster_ageRange") := .(NA, NA, NA, NA, NA)]
dat_cases[, status := "UNV"]
dat_cases[date > date_1, status := "PAR"]
dat_cases[date > vax_date, status := "VAX"]
dat_cases[date > booster_date, status := "BST"]
dat_cases[date <= date_1, manu := "UNV"]
dat_cases$status <- factor(dat_cases$status, levels = c("UNV", "PAR", "VAX", "BST"))
dat_cases <- dat_cases[date <= analysis_last_day,]

unvax <- dat_cases[status == "UNV" & gender %in% c("F", "M") & !ageRange %in% c("0-11", "12-17") &
                     date <= analysis_last_day, ]
unvax <- unvax[, .(cases = .N, hosp = sum(hosp), death = sum(death)), 
               keyby = .(date, ageRange, gender)]
unvax <- merge(pop_unvax, unvax, by = c("date", "ageRange", "gender"), all.x = TRUE) 
unvax[is.na(unvax)] <- 0
unvax <- melt(unvax, measure.vars = c("cases", "hosp", "death"),
              variable.name = "outcome", value.name = "obs")
unvax[, outcome := factor(outcome, levels = c("cases", "hosp","death"))]
unvax[, ageRange := as.character(ageRange)]
unvax[outcome == "death", 
      ageRange := as.character(fct_collapse(
        factor(ageRange), "18-44" = c("18-24","25-34","35-44")))]
unvax = unvax[,.(obs=sum(obs), poblacion = sum(poblacion)), 
              keyby = c("date","ageRange","gender", "outcome")]   

vax <- dat_cases[status == "VAX" & gender %in% c("F", "M") & !ageRange %in% c("0-11", "12-17") &
                   date <= analysis_last_day, ]
vax[,day := as.numeric(date) - as.numeric(vax_date)]
vax <- vax[day <= the_ndays, ]
vax <- vax[, .(cases = .N, hosp = sum(hosp), death = sum(death)), 
             keyby = .(manu, ageRange, gender, date, vax_date)]
vax <- merge(pop_vax, vax, by = c("manu", "ageRange", "gender", "date", "vax_date"), all.x = TRUE) 
vax[is.na(vax)] <- 0
vax <- melt(vax, measure.vars = c("cases", "hosp", "death"),
            variable.name = "outcome", value.name = "obs")
vax[, outcome := factor(outcome, levels = c("cases", "hosp","death"))]
vax[, ageRange := as.character(ageRange)]
vax[outcome == "death", 
    ageRange := as.character(fct_collapse(
      factor(ageRange), "18-44" = c("18-24","25-34","35-44")))]
vax[,day:=as.numeric(date-vax_date)]
vax <- vax[,.(obs=sum(obs), poblacion = sum(poblacion), day = day[1]), 
  keyby = c("manu", "ageRange", "gender", "date", "vax_date", "outcome")]


### booster counts by date
all_dates <-  seq(first_booster_day, analysis_last_day, "days")
all_combs <- CJ(date = all_dates,
                booster_date = all_dates, 
                ageRange = levels(dat_vax$ageRange_1)[-c(1,2)],
                gender = c("M", "F"),
                manu_1 = levels(dat_vax$manu_1),
                booster_manu = levels(dat_vax$booster_manu))
all_combs <- all_combs[as.numeric(date - booster_date) >=0 & 
                         as.numeric(date - booster_date) <= the_booster_ndays,]

combo_map <- CJ(manu = manu_levels[-1],
                booster_manu = manu_levels[-1])
combo_map$manu <- factor(combo_map$manu, manu_levels)
combo_map$booster_manu <- factor(combo_map$booster_manu, manu_levels)
combo_map <- combo_map[order(manu, booster_manu)]
combo_map$combo <- with(combo_map, 
                        case_when(manu == booster_manu &  manu != "JSN" ~ paste(manu,"+",booster_manu),
                                  manu != "JSN" & booster_manu !="JSN" ~ "MOD + PFR or PFR + MOD",
                                  manu == "JSN" & booster_manu !="JSN" ~ "JSN + MOD or PFR",
                                  TRUE ~ "JSN booster"))

  
bst <- dat_cases[status == "BST" & gender %in% c("F", "M") & !ageRange %in% c("0-11", "12-17") &
                   date <= analysis_last_day, ]
bst <- bst[, .(cases = .N, hosp = sum(hosp), death = sum(death)), 
           keyby = .(manu_1, booster_manu, ageRange, gender, date, booster_date)]
bst <- merge(all_combs, bst, 
             by = c("date","booster_date", "ageRange", "gender", "manu_1", "booster_manu"), 
             all.x = TRUE) 
bst[is.na(bst)] <- 0
bst <- merge(bst, daily_counts_booster_age_gender_manu,
             by = c("manu_1", "booster_manu", "ageRange", "gender", "booster_date"), 
             all.x = TRUE) 
bst <- bst[!is.na(poblacion),]
setnames(bst, c("manu_1", "booster_date"), c("manu", "vax_date"))
bst <- bst[, manu := factor(manu, levels = manu_levels)]
bst <- melt(bst, measure.vars = c("cases", "hosp", "death"),
            variable.name = "outcome", value.name = "obs")
bst[, outcome := factor(outcome, levels = c("cases", "hosp","death"))]
bst[, ageRange := as.character(ageRange)]
bst[outcome == "death", 
    ageRange := as.character(fct_collapse(
      factor(ageRange), "18-44" = c("18-24","25-34","35-44")))]
bst[,day:=as.numeric(date-vax_date)]
bst <- bst[,.(obs=sum(obs), poblacion = sum(poblacion), day = day[1]), 
           keyby = c("manu", "booster_manu","ageRange", "gender", "date", "vax_date", "outcome")]
bst <- merge(bst, combo_map, by = c("manu", "booster_manu"), all.x = TRUE)


source("funcs.R")
library("splines")

exp <- unvax[, compute_expected(.SD), keyby = "outcome"]
exp[, rate := fit / poblacion]
exp <- exp[,c("outcome", "date", "ageRange", "gender", "rate")]

#exp %>% ggplot(aes(date, rate, color = gender)) + geom_line() +facet_grid(outcome~ageRange, scales="free_y")

vax_waning_fits <- vax %>%
  left_join(exp, by = c("outcome", "date", "ageRange", "gender")) %>%
  filter(poblacion>0) %>%
  mutate(exp = rate * poblacion,
         day = ifelse(manu=="JSN", pmin(day, 6*30), pmin(day, 7*30))) %>%
  group_by(outcome, manu, day) %>%
  summarize(obs = sum(obs), exp = sum(exp), .groups = "drop") %>%
  group_by(outcome, manu) %>%
  do(fit_waning_model(., max.day = 240, knots=c(60)))


bst_waning_fits <- bst %>%
  filter(outcome == "cases") %>%
  left_join(exp, by = c("outcome", "date", "ageRange", "gender")) %>%
  filter(poblacion>0) %>%
  mutate(exp = rate * poblacion,
         day = ifelse(manu=="JSN", pmin(day, 2*30), pmin(day, 3*30))) %>%
  group_by(manu, day) %>%
  summarize(obs = sum(obs), exp = sum(exp), .groups = "drop") %>%
  group_by(manu) %>%
  do(fit_waning_model_linear(., max.day = 90))



## make tables
tmp <- vax %>%
  left_join(exp, by = c("outcome", "date", "ageRange", "gender")) %>%
  filter(poblacion>0) %>%
  mutate(exp = rate * poblacion)

vax_cases_tab  <- tmp %>%
  filter(outcome=="cases") %>%
  mutate(exp = poblacion*rate, months = pmin(floor(day/30/2)*2+2, 8)) %>%
  group_by(outcome, manu, months) %>%
  summarize(obs=sum(obs), exp=sum(exp), .groups = "drop") %>%
  mutate(eff = 1-obs/exp, conf.low = 1 - qpois(0.975, obs, )/exp, 
         conf.high = 1-qpois(0.025, obs,)/exp)  %>%
  bind_rows(vax_waning_fits %>% 
              filter(outcome=="cases" & day==0) %>% 
              rename(eff = fit, months = day)) %>%
  arrange(months, manu) %>%
  mutate(months = recode(as.character(months), "2"="(0-2]", "4"="(2-4]", "6"="(4-6]","8" = "6+")) %>%
  mutate(text = paste0(make_pct(eff), " (", make_pct(conf.low), "-",make_pct(conf.high), ")")) %>%
  select(manu, months, text)

## Events
tmp_1 <- tmp %>% 
  filter(manu != "JSN" & outcome!="cases") %>%
  mutate(exp = poblacion*rate, months = ifelse(day/30<6, "(0-6]", "6+")) 
tmp_2 <- tmp %>% 
  filter(manu == "JSN"& outcome!="cases") %>%
  mutate(exp = poblacion*rate, months = ifelse(day/30<2, "(0-2]", "2+")) 

vax_events_tab <- bind_rows(tmp_1, tmp_2) %>%
  group_by(outcome, manu, months) %>%
  summarize(obs=sum(obs), exp=sum(exp), .groups = "drop") %>%
  mutate(eff = 1-obs/exp, conf.low = 1 - qpois(0.975, obs, )/exp, 
         conf.high = 1-qpois(0.025, obs,)/exp)  %>%
  bind_rows(vax_waning_fits %>% 
              filter(outcome!="cases" & day==0) %>% 
              mutate(day=as.character(day))%>%
              rename(eff = fit, months = day)) %>%
  mutate(text = paste0(make_pct(eff), " (", make_pct(conf.low), "-",make_pct(conf.high), ")")) %>%
  mutate(months = factor(months, levels = c("0", "(0-6]","(0-2]","6+","2+")))%>%
  arrange(outcome, months, manu) %>%
  select(outcome, manu, months, text)


booster_totals <- copy(daily_counts_booster_age_gender_manu)
booster_totals <- setnames(booster_totals, "manu_1", "manu")
booster_totals <- merge(booster_totals, combo_map, by =  c("manu", "booster_manu"), all.x = TRUE)
booster_full_totals <- booster_totals[ , .(total = sum(poblacion)), keyby = c("manu", "booster_manu")]
booster_totals <- booster_totals[ , .(total = sum(poblacion)), keyby = "combo"]

bst_cases_tab <- bst %>% filter(outcome == "cases") %>%
  left_join(exp, by = c("outcome", "date", "ageRange", "gender")) %>%
  filter(poblacion>0) %>%
  mutate(exp = poblacion*rate) %>%
  group_by(combo) %>%
  summarize(obs=sum(obs), exp=sum(exp), .groups = "drop") %>%
  mutate(eff = 1-obs/exp, conf.low = 1 - qpois(0.975, obs, )/exp, conf.high = 1-qpois(0.025, obs,)/exp)  %>%
  mutate(text = paste0(make_pct(eff), " (", make_pct(conf.low), "-",make_pct(conf.high), ")"))%>%
  filter(combo != "JSN booster") %>%
  select(combo, text)

bst_events_tab <- bst %>% filter(outcome != "cases") %>%
  left_join(exp, by = c("outcome", "date", "ageRange", "gender")) %>%
  filter(poblacion>0) %>%
  mutate(exp = poblacion*rate) %>%
  group_by(outcome) %>%
  summarize(obs=sum(obs), exp=sum(exp), .groups = "drop") %>%
  mutate(eff = 1-obs/exp, conf.low = 1 - qpois(0.975, obs, )/exp, conf.high = 1-qpois(0.025, obs,)/exp)  %>%
  mutate(text = paste0(make_pct(eff), " (", make_pct(conf.low), "-",make_pct(conf.high), ")")) %>%
  select(outcome, text)

## Not enough data yet
if(FALSE){
cases <- vax[outcome == "cases", ]
cases <- cases[, after_booster_date := ifelse(manu=="JSN", day>=2*30, day>=6*30)]
cases <- cases[,.(poblacion = sum(obs)), keyby =.(manu,  ageRange, after_booster_date)]

events <- vax[outcome != "cases", !"poblacion"]
events <- events[,after_booster_date := ifelse(manu=="JSN", day>=2*30, day>=6*30)]
events <- events[,.(obs = sum(obs)), keyby =.(manu,  ageRange, after_booster_date, outcome)]
events <- merge(cases, events, by = c("manu",  "ageRange", "after_booster_date"), all.x=TRUE)
events[is.na(events)] <-0
events[,rate := obs/poblacion]


cases <- bst[outcome == "cases", ]
cases <- cases[,.(poblacion = sum(obs)), keyby =.(manu,  ageRange)]

bst_events <- bst[outcome != "cases", !"poblacion"]
bst_events <- bst_events[,.(obs = sum(obs)), keyby =.(manu,  ageRange, outcome)]
bst_events <- merge(cases, bst_events, by = c("manu",  "ageRange"), all.x=TRUE)
bst_events[is.na(bst_events)] <-0
bst_events[,rate := obs/poblacion]
}

totals <- list(vaxed = sum(!is.na(dat_vax$vax_date)),
                 boosted = sum(!is.na(dat_vax$booster_date)),
               cases = sum(dat_cases_vax$date>=first_booster_day),
               all_cases = nrow(dat_cases))
               
save(pr_pop, first_booster_day, first_day, analysis_last_day, totals, vax_waning_fits, bst_events_tab, bst_cases_tab, booster_full_totals, bst_waning_fits,
     booster_totals, vax_cases_tab, vax_events_tab,  file = "rdas/paper-data.rda")

