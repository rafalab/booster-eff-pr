make_pct <- function(x, digit = 0) ifelse(is.na(x), "", paste0(format(round(100*x, digit = digit), nsmall = digit), "%"))

make_pretty <- function(x) prettyNum(replace_na(x, " "), big.mark = ",")

dynamic_round <- function(x, min_rounded = 10){
  ifelse(round(x) >= min_rounded, 
         prettyNum(round(x), big.mark = ","),
         prettyNum(round(x, digits = 1), nsmall = 1, big.mark = ",")) %>%
    replace_na("")
}

compute_expected <- function(tab, days.between.knots = 45){
  df <- round(length(unique(tab$date))/days.between.knots)
  
  dat <- tab %>%
    mutate(x = as.numeric(date) - as.numeric(min(date)),
           wd = factor(wday(date)))
  
  fit <- glm(obs ~ gender + wd + ageRange*ns(x, df = df), offset = log(poblacion), 
             family = "poisson", data = dat)
  
  pred <- predict(fit, se.fit = TRUE, type = "response")
  tab <- mutate(tab, fit = pred$fit, se.fit = pred$se.fit)
  return(tab)  
}

fit_waning_model <- function(tab, alpha = 0.05, max.day = 240, knots = 120, transf = function(x){ 1-(1/exp(-x)) }){
  
  tab <- filter(tab, day <= max.day)
  
  fit <- glm(obs ~ ns(day, knots = knots), offset = log(exp), family = "quasipoisson", data = tab)
  
  pred <- predict(fit, newdata = list(day=tab$day, exp = 1), se = TRUE)
  
  mutate(tab, 
         fit = transf(pred$fit),
         conf.low = transf(pred$fit + qnorm(1-alpha/2) * pred$se.fit),
         conf.high = transf(pred$fit - qnorm(1-alpha/2) * pred$se.fit))
}

fit_waning_model_linear <- function(tab, alpha = 0.05, max.day = 240,  transf = function(x){ 1-(1/exp(-x)) }){
  
  tab <- filter(tab, day <= max.day)
  
  fit <- glm(obs ~ day, offset = log(exp), family = "quasipoisson", data = tab)
  
  pred <- predict(fit, newdata = list(day=tab$day, exp = 1), se = TRUE)
  
  mutate(tab, 
         fit = transf(pred$fit),
         conf.low = transf(pred$fit + qnorm(1-alpha/2) * pred$se.fit),
         conf.high = transf(pred$fit - qnorm(1-alpha/2) * pred$se.fit))
}
