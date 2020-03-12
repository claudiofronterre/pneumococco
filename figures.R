# Load packages ----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("ggplot2", "dplyr", "ggsci", "JLutils", "extrafont", "MASS") # package names
pacman::p_load(pkgs, character.only = T)

# Model output -----------------------------------------
source("R/censorModel2.R")
cocco <- readRDS("data/cocco_clean.rds")
cocco <- cocco[!(cocco$pcv_vaccd == "PCV13-unvacc'd" & cocco$surv == 4 & cocco$age_y < 1), ] 
tcu <- min(cocco$age_y[cocco$age_y >= 0.5 & cocco$pcv_vaccd == 0])
tcv <- tcu
fit <- censorModel(formula = vt_01 ~ 1, 
                   data = cocco, 
                   ini.par = c(0.2, 0.05, 3), 
                   vacc = "pcv_vaccd",
                   age = "age_y",
                   tcu = tcv,
                   tcv = tcv,
                   method = "Nelder-Mead",
                   epsilon = 0.0001)

library(MASS)
mu <- as.numeric(fit$par)
covar <- fit$covar
theta <- mvrnorm(n = 1000, mu = mu, Sigma = covar)

logistic <- function(x) 1 / (1 + exp(- x))
p <- ncol(theta) - 3
thetau <- cbind(theta[, 1:p], exp(theta[, p + 1]), theta[, p + 3])
thetav <- cbind(theta[, 1:p], exp(theta[, p + 2]))

D <- fit$D
age <- seq(0, 11, l = 1000)

predictfu <- function(x, age) {
  xi <- x[1:p]
  alpha <- logistic(D %*% xi)
  alpha <- mean(alpha)
  delta <- x[p + 1]
  beta <- x[p + 2]
  alpha * beta * exp(- delta * (age - tcv))
} 

predictfv <- function(x, age) {
  xi <- x[1:p]
  alpha <- logistic(D %*% xi)
  alpha <- mean(alpha)
  delta <- x[p + 1]
  alpha * exp(- delta * (age - tcv))
} 


fit_unvacc <- t(apply(apply(thetau, 1, predictfu, age = age), 1, 
                      quantile, prob = c(0.025, 0.975)))
fit_vacc <- t(apply(apply(thetav, 1, predictfv, age = age), 1, 
                    quantile, prob = c(0.025, 0.975)))

fit.vacc.low <- fit_vacc[, 1]
fit.vacc <- apply(apply(thetav, 1, predictfv, age = age), 1, mean)
fit.vacc.up <- fit_vacc[, 2]

fit.unvacc.low <- fit_unvacc[, 1]
fit.unvacc <- apply(apply(thetau, 1, predictfu, age = age), 1, mean)
fit.unvacc.up <- fit_unvacc[, 2]

df <- data.frame(cat = as.factor(rep(c("vaccinated", "unvaccinated"), each = length(age))),
                 probability = c(fit.vacc, fit.unvacc),
                 low = c(fit.vacc.low, fit.unvacc.low),
                 up = c(fit.vacc.up, fit.unvacc.up),
                 age = c(age, age))

pal <- get_palette(palette = "lancet", 2)
size <- 11
legend <- "top"
xlim <- 0:11

df[df$age < tcu & df$cat == "unvaccinated", -1] <- NA 
df[df$age < tcv & df$cat == "vaccinated", -1] <- NA 

ggplot(data = df, aes(x = age, y = probability, color = cat)) + 
  geom_line(size = 0.8) +
  xlab("Age (years)") +
  ylab("Probability of VT carriage") +
  scale_color_manual(values = rev(pal), name = "", lab = c("Unvaccinated", "Vaccinated")) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = cat), alpha = 0.2, show.legend = F,
              color = NA) +
  theme_bw(base_size = size) +
  theme(legend.position = legend, text = element_text(family = "Arial", size = 12)) +
  scale_x_continuous(breaks = xlim, limits = c(0, 11))

ggsave("figs/model.tiff", dpi = 300, units = "cm", width = 20, height = 14)

# Age distribution of children by survey -----------------
cocco <- readr::read_csv("data/PCVPA_31dec_adj_claudio_.csv")
cocco <- cocco[!(cocco$pcv_vaccd == "PCV13-unvacc'd" & cocco$surv == 4 & cocco$age_y < 1), ]
cocco <- cocco[-which(is.na(cocco$surv)), ] # remove NA


cocco %>% 
  filter(typ_cat7 != "adult") %>%
  ggplot(aes(x = age_y, col = pcv_vaccd, fill = pcv_vaccd)) +
  geom_histogram(position = "identity", alpha = .5, show.legend = T) +
  facet_wrap(~ surv, labeller = labeller(surv = function(x) paste("Survey", x))) +
  scale_x_continuous(breaks = round(seq(0, 11, by = 1))) +
  scale_color_manual(values = rev(pal_lancet()(2)), 
                     name = "", label = c("Unvaccinated", "Vaccinated")) +
  scale_fill_manual(values = rev(pal_lancet()(2)), 
                    name =  "", label = c("Unvaccinated", "Vaccinated")) +
  labs(y = "Count", 
       x = "Age (years)") +
  theme_bw() +
  theme(legend.position = c(.89, .17), text = element_text(family = "Arial", size = 12))
ggsave("figs/age_survey.tiff", dpi = 300, units = "cm", width = 20, height = 14)

# Ditribution of pneumococco carriage over different age groups -----------------
cocco$age_cat <- cut(cocco$age_y, breaks = seq(0, 11, by = 1))
cocco$pcv_vaccd <- ifelse(cocco$pcv_vaccd == "PCV13-vacc'd", "Vaccinated", "Unvaccinated")
cocco %>% 
  filter(typ_cat7 != "adult") %>%
  ggplot(aes(x = age_cat, fill = factor(vt_01))) +
  geom_bar(position = "fill", show.legend = T) +
  stat_fill_labels(size = 3) +
  facet_wrap(~ pcv_vaccd, ncol = 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(name = "Vaccine type carriage", labels = c("No", "Yes"),
                    type = "qual", palette = 7,
                    position = "top") +
  labs(x = "Age (years)", y = "") +
  theme_bw() +
  theme(legend.position = "top", 
        panel.grid = element_blank(), 
        text = element_text(family = "Arial", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/coccodist.tiff", dpi = 300, units = "cm", width = 20, height = 14)


# Residuals ----------------------------------------
source("R/plotResiduals.R")
plotResiduals(model = fit, ngroups = 15, legend = "top") +
  coord_equal() +
  labs(y = "Residuals", x = "Predicted") +
  theme(text = element_text(family = "Arial", size = 12),
        panel.grid = element_blank()) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, col = "grey")
ggsave("figs/residuals.tiff", dpi = 300, units = "cm", width = 20, height = 7)
