# Load required packages and functions
if (!require("pacman")) install.packages("pacman")
pkgs = c("dplyr") # package names
pacman::p_load(pkgs, character.only = T)
source("R/coccoModel.R")

# Load data
cocco <- readRDS("data/cocco_clean.rds")

# Model with same delta (same rate of decay for vacc and unvacc)
fit <- coccoModel(response = "vt_01", vacc = "pcv_vaccd", age = "age_y",
                  data = cocco, delta2 = F, 
                  kv = 0.3452055, ku = 1,
                  ini.par = c(0.5, 0.1))
fit$results
fit$lik

fit2 <- coccoModel(response = "vt_01", vacc = "pcv_vaccd", age = "age_y",
                   data = cocco, delta2 = T, 
                   kv = 0.3452055, ku = 1,
                   ini.par = c(-1, 3, 0.1, 0.1))
fit2$results
fit2$lik

age0 <- 1 / 12
age1 <- (1 / 12) * 1
beta <- 2.66
abs <- (psych::logistic(beta * age1) - 0.5) - (psych::logistic(beta * age0) - 0.5)
rel <- abs / (psych::logistic(beta * age0) - 0.5)


source("R/coccoModel2.R")
fit3 <- coccoModel2(response = "vt_01", vacc = "pcv_vaccd", age = "age_y",
                    data = cocco, delta2 = T, 
                    kv = 0.3452055, ku = 1,
                    ini.par = c(0.5, -1, 1, 0.1, 0.1))
fit3$results
fit3$lik



# Plot the fitted model
source("R/plot_model.R")
plot_model(fit, size = 11, xlim = 0:11, legend = "top") +
  xlab("age (years)") +
  

source("R/plot_model2.R")
plot_model2(fit2, size = 11, age = seq(0, 10, l = 1000), xlim = 0:10, legend = "top") +
  xlab("age (years)") +
  labs(title = "Natural immunity set point = 24 months")

source("R/plot_model3.R")
plot_model3(fit3, size = 11, age = seq(0, 10, l = 1000), xlim = 0:10, legend = "top") +
  xlab("age (years)")
# Plot residuals

source("R/plotResiduals.R")

plotResiduals(model = fit, k = 90, size = 11, legend = "top")


# Bootstrap to obtain cofidence bands for predicted curve
library(MASS)
mu <- as.numeric(fit2$par)
covar <- fit2$covar
theta <- mvrnorm(n = 1000, mu = mu, Sigma = covar)

thetau <- cbind(theta[, 1:2], exp(theta[, 3]), fit2$ku)
thetav <- cbind(theta[, 1:2], exp(theta[, 4]), fit2$kv)

predictf <- function(x, age) {
  alpha <- x[1]
  beta <- x[2]
  delta <- x[3]
  k <- x[4]
  logistic <- function(x) 1 / (1 + exp(- x))
  logistic(alpha + beta * age) * (age < k) + 
    logistic(alpha + beta * k) * exp(- delta * (age - k)) * (age >= k)
} 

age <- seq(0, 10, l = 1000)

fit_unvacc <- t(apply(apply(thetau, 1, predictf, age = age), 1, 
                      quantile, prob = c(0.025, 0.975)))
fit_vacc <- t(apply(apply(thetav, 1, predictf, age = age), 1, 
                    quantile, prob = c(0.025, 0.975)))

mean_par_vacc <- apply(thetav, 2, mean)
mean_par_unvacc <- apply(thetau, 2, mean)

fit.vacc <- predictf(mean_par_vacc, age)
fit.vacc.low <- fit_vacc[, 1]
fit.vacc.up <- fit_vacc[, 3]

fit.unvacc <- predictf(mean_par_unvacc, age)
fit.unvacc.low <- fit_unvacc[, 1]
fit.unvacc.up <- fit_unvacc[, 3]


df <- data.frame(cat = as.factor(rep(c("vaccinated", "unvaccinated"), each = length(age))),
                 probability = c(fit.vacc, fit.unvacc),
                 low = c(fit.vacc.low, fit.unvacc.low),
                 up = c(fit.vacc.up, fit.unvacc.up),
                 age = c(age, age))

require(ggplot2)
require(ggpubr)

pal <- get_palette(palette = "lancet", 2)
size <- 11
legend <- "top"
xlim <- 0:10

ggplot(data = df, aes(x = age, y = probability, color = cat)) + 
  geom_line(size = 1.2) +
  xlab("age") +
  ylab("probability of VT carriage") +
  scale_color_manual(values = rev(pal), name = "") +
  geom_ribbon(aes(ymin = low, ymax = up, fill = cat), alpha = 0.2, show.legend = F,
              color = NA) +
  theme_bw(base_size = size) +
  theme(legend.position = legend) +
  scale_x_continuous(breaks = xlim) 


# Profile likelihood
plik_beta <- function(beta, ini.par = c(0.1, 1, 2)) {
  par.start <- log(ini.par)
  log.lik <- function(par) {
    # Model parameters
    delta <- exp(par[1])
    kv <- exp(par[2])
    ku <- exp(par[3])
    
    # Probability functional form (vaccinated)
    pr[vacc] <- beta * age[vacc] * (age[vacc] < kv) + 
      (beta * kv) * exp(- delta * (age[vacc] - kv)) * (age[vacc] >= kv)
    
    # Probability functional form (vaccinated)
    pr[unvacc] <- beta * age[unvacc] * (age[unvacc] < ku) + 
      (beta * ku) * exp(- delta * (age[unvacc] - ku)) * (age[unvacc] >= ku)
    
    # Negative log-likelihood
    if(any(pr <= 0 | pr >= 1)) {
      out <- 1500000
    } else {
      out <- - (sum(log(1 - pr)) + sum(y * log(pr / (1 - pr))))
    }    
    as.numeric(out)
  }
  
  estim1 <- list(); estim1$value <- 1
  estim2 <- list(); estim2$value <- 2
  
  epsilon <- 0.00000001
  method <- "Nelder-Mead"
  
  while(abs(estim2$value - estim1$value) > epsilon) {
    estim1 <- optim(par = par.start, fn = log.lik, hessian = T, control = list(maxit = 5000), method = method)  
    estim2 <- optim(par = estim1$par, fn = log.lik, hessian = T, control = list(maxit = 5000), method = method)
    par.start <- estim2$par
  }
  - estim2$value
}

beta <- seq(0, 1, l = 10)
plik <- sapply(beta, function(x) plik_beta(x))  

plot(beta, plik, type = "l")
t <- .5
plot(beta[beta < t], plik[beta < t], type = "l")

beta[which.max(plik)]

exp(plik_beta(beta[which.max(plik)]))
