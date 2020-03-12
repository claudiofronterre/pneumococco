library(brms)
theme_set(theme_default())

library(psych)

# Simulate data

n <- 2000
delta <- 0.1
beta <- 0.5
alpha <- 0.7
gamma <- 0.2

v <- sample(c(1, 0), n, replace = T)
age <- rep(0, n) ; t <- rep(0, n)
age[v == 1] <- runif(sum(v == 1), 0, 6)
age[v == 0] <- runif(sum(v!= 1), 5, 10)
t[v == 1] <- runif(sum(v == 1), 0.1, 2)

p <- alpha * beta ^ v * exp(- delta * age - v * gamma * t)
y <- rbinom(n, size = 1, prob = p)

f <- bf(y ~ alpha * beta ^ v * exp(- delta * age - v * gamma * t), 
        alpha + beta + delta + gamma ~ 1, 
        nl = TRUE)

priors <- prior(beta(1, 1), nlpar = "alpha", lb = 0, ub = 1) +
  prior(beta(1, 1), nlpar = "beta", lb = 0, ub = 1) +
  prior(cauchy(0, 10), nlpar = "delta", lb = 0) +
  prior(cauchy(0, 10), nlpar = "gamma", lb = 0) 

df <- data.frame(y, v, age, t)
df$v <- as.factor(df$v)

fit1 <- decayModel(formula = y ~ 1, 
                   data = df, 
                   ini.par = c(0.5, 0.1, 0.3), 
                   vacc = "v",
                   age = "age",
                   agev = "t",
                   tc = 0,
                   gamma = F,
                   delta2 = F, 
                   epsilon = 0.00000001, 
                   eig = F)

fit <- brm(formula = f, family = bernoulli("identity"), prior = priors, 
           data = df, cores = 4)
summary(fit)
plot(fit)
gg <- marginal_effects(fit, effects = "age:v")
plot(gg, plot = F)[[1]] +
  ggsci::scale_colour_lancet(name = "Vaccinated", labels = c("no", "yes"), position = "top") +
  ggsci::scale_fill_lancet(name = "Vaccinated", labels = c("no", "yes"), position = "top")
