plot_model3 <- function(model, size, age, xlim, legend = "bottom") {
  
  require(ggplot2)
  require(ggpubr)
  
  beta <- model$res["beta", 1]; beta.low <- model$res["beta",2]; beta.up <- model$res["beta",3]
  deltau <- model$res["deltau", 1]; deltau.low <- model$res["deltau",2]; deltau.up <- model$res["deltau",3]
  deltav <- model$res["deltav", 1]; deltav.low <- model$res["deltav",2]; deltav.up <- model$res["deltav",3]
  gamma1 <- model$res["gamma1", 1]; gamma1.low <- model$res["gamma1",2]; gamma1.up <- model$res["gamma1",3]
  gamma2 <- model$res["gamma2", 1]; gamma2.low <- model$res["gamma2",2]; gamma2.up <- model$res["gamma2",3]
  kv <- model$kv
  ku <- model$ku
  logistic <- function(x) 1 / (1 + exp(- x))
  
  fit.unvacc.low <- (logistic(beta.low * age) - 0.5) * (age < ku) + 
    (logistic(beta.low * ku) - 0.5) * exp(- deltau.up * (age - ku)) * (age >= ku)
  fit.unvacc <- (logistic(beta * age) - 0.5) * (age < ku) + 
    (logistic(beta * ku) - 0.5) * exp(- deltau * (age - ku)) * (age >= ku)
  fit.unvacc.up <- (logistic(beta.up * age) - 0.5) * (age < ku) + 
    (logistic(beta.up * ku) - 0.5) * exp(- deltau.low * (age - ku)) * (age >= ku)
  
  fit.vacc.low <- (logistic(beta.low * age) - 0.5) * (age < kv) + 
    (logistic(beta.low * kv) - 0.5) * (0.5 + logistic(gamma1.low * (age - kv) + gamma2.low * (age - kv) ^ 2)) * (age >= kv & age < ku) + (logistic(beta.low * kv) - 0.5) * (0.5 + logistic(gamma1.low * (ku - kv) + gamma2.low * ((ku - kv) ^ 2))) * exp(- deltav.up * (age - ku)) * (age >= ku)
  fit.vacc <- (logistic(beta * age) - 0.5) * (age < kv) + 
    (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (age - kv) + gamma2 * (age - kv) ^ 2)) * (age >= kv & age < ku) + (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (ku - kv) + gamma2 * ((ku - kv) ^ 2))) * exp(- deltav * (age - ku)) * (age >= ku)
  fit.vacc.up <- (logistic(beta.up * age) - 0.5) * (age < kv) + 
    (logistic(beta.up * kv) - 0.5) * (0.5 + logistic(gamma1.up * (age - kv) + gamma2.up * (age - kv) ^ 2)) * (age >= kv & age < ku) + (logistic(beta.up * kv) - 0.5) * (0.5 + logistic(gamma1.up * (ku - kv) + gamma2.up * ((ku - kv) ^ 2))) * exp(- deltav.low * (age - ku)) * (age >= ku)
  
  
  df <- data.frame(cat = as.factor(rep(c("vaccinated", "unvaccinated"), each = length(age))),
                   probability = c(fit.vacc, fit.unvacc),
                   low = c(fit.vacc.low, fit.unvacc.low),
                   up = c(fit.vacc.up, fit.unvacc.up),
                   age = c(age, age))
  
  pal <- get_palette(palette = "lancet", 2)
  
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
    
}
