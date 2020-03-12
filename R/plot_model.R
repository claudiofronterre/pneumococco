plot_model <- function(model, size, xlim, legend = "bottom") {
  
  require(ggplot2)
  require(ggpubr)
  
  delta2 <- nrow(model$results) > 2
  
  if (delta2) {
    
  } else {
    beta <- model$res["beta", 1]; beta.low <- model$res["beta",2]; beta.up <- model$res["beta",3]
    delta <- model$res["delta", 1]; delta.low <- model$res["delta",2]; delta.up <- model$res["delta",3]
    kv <- model$kv
    ku <- model$ku
    age <- seq(0, 11, l = 1000)
    logistic <- function(x) 1 / (1 + exp(- x))
    
    fit.unvacc.low <- (logistic(beta.low * age) - 0.5) * (age < ku) + 
      (logistic(beta.low * ku) - 0.5) * exp(- delta.up * (age - ku)) * (age >= ku)
    fit.unvacc <- (logistic(beta * age) - 0.5) * (age < ku) + 
      (logistic(beta * ku) - 0.5) * exp(- delta * (age - ku)) * (age >= ku)
    fit.unvacc.up <- (logistic(beta.up * age) - 0.5) * (age < ku) + 
      (logistic(beta.up * ku) - 0.5) * exp(- delta.low * (age - ku)) * (age >= ku)
    
    fit.vacc.low <- (logistic(beta.low * age) - 0.5) * (age < kv) + 
      (logistic(beta.low * kv) - 0.5) * exp(- delta.up * (age - kv)) * (age >= kv)
    fit.vacc <- (logistic(beta * age) - 0.5) * (age < kv) + 
      (logistic(beta * kv) - 0.5) * exp(- delta * (age - kv)) * (age >= kv)
    fit.vacc.up <- (logistic(beta.up * age) - 0.5) * (age < kv) + 
      (logistic(beta.up * kv) - 0.5) * exp(- delta.low * (age - kv)) * (age >= kv)
  }
  
  
  
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
