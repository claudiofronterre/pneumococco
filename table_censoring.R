source("R/censorModel2.R")
cocco <- readRDS("data/cocco_clean.rds")
cocco <- cocco[!(cocco$pcv_vaccd == "PCV13-unvacc'd" & cocco$surv == 4 & cocco$age_y < 1), ] 
tcu <- min(cocco$age_y[cocco$age_y >= 0.5 & cocco$pcv_vaccd == 0])
tcv <- seq(4, 7, .5)

for (i in 1:length(tcv)) {
  fit <- censorModel(formula = vt_01 ~ 1, 
                     data = cocco, 
                     ini.par = c(0.2, 0.05, 3), 
                     vacc = "pcv_vaccd",
                     age = "age_y",
                     tcu = tcv[i],
                     tcv = tcv[i],
                     method = "Nelder-Mead",
                     epsilon = 0.0001)
  print(tcv[i])
  print(fit$results)
}
