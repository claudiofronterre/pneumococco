coccoModel2 <- function(response,
                       vacc,
                       age,
                       data, 
                       ini.par, 
                       ku,
                       kv,
                       delta2,
                       method = "Nelder-Mead",
                       epsilon = 0.0001) {
  y <- cocco[[response]]
  vaccinated_dummy <- cocco[[vacc]]
  vacc <- vaccinated_dummy == 1
  unvacc <- vaccinated_dummy == 0
  age <- cocco[[age]]
  pr <- rep(0, length(y)) 
  logistic <- function(x) 1 / (1 + exp(- x))
  
  par.start <- c(ini.par[1:3], log(ini.par[4:length(ini.par)])) 
  
  if(delta2) {
    log.lik <- function(par) {
      # Model parameters
      beta <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      deltau <- exp(par[4])
      deltav <- exp(par[5])
      
      # Probability functional form (vaccinated)
      pr[vacc] <- (logistic(beta * age[vacc]) - 0.5) * (age[vacc] < kv) + 
        (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (age[vacc] - kv) + gamma2 * (age[vacc] - kv) ^ 2)) * (age[vacc] >= kv & age[vacc] < ku) + (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (ku - kv) + gamma2 * ((ku - kv) ^ 2))) * exp(- deltav * (age[vacc] - ku)) * (age[vacc] >= ku)
      
      # Probability functional form (vaccinated)
      pr[unvacc] <- (logistic(beta * age[unvacc]) - 0.5) * (age[unvacc] < ku) + 
        (logistic(beta * ku) - 0.5) * exp(- deltau * (age[unvacc] - ku)) * (age[unvacc] >= ku)
      
      # Negative log-likelihood
      out <- - (sum(log(1 - pr)) + sum(y * log(pr / (1 - pr))))
      as.numeric(out)
    }
    } else {
      log.lik <- function(par) {
        # Model parameters
        beta <- par[1]
        delta <- exp(par[2])
        
        # Probability functional form (vaccinated)
        pr[vacc] <- (logistic(beta * age[vacc]) - 0.5) * (age[vacc] < kv) + 
          (logistic(beta * kv) - 0.5) * exp(- delta * (age[vacc] - kv)) * (age[vacc] >= kv)
        
        # Probability functional form (vaccinated)
        pr[unvacc] <- (logistic(beta * age[unvacc]) - 0.5) * (age[unvacc] < ku) + 
          (logistic(beta * ku) - 0.5) * exp(- delta * (age[unvacc] - ku)) * (age[unvacc] >= ku)
        
        # Negative log-likelihood
        out <- - (sum(log(1 - pr)) + sum(y * log(pr / (1 - pr))))
        as.numeric(out)
    }
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
  
  
  par <- c(estim2$par[1:3], exp(estim2$par[4:length(ini.par)]))
  if (delta2) {
    names(par) <- c("beta", "gamma1", "gamma2", "deltau", "deltav")
  } else {
    names(par) <- c("beta", "delta")
  }
  
  ## Confidence intervals
  conf.level <- 0.95
  z <- - qnorm((1 - conf.level) / 2)
  
  conf.int <- cbind(estim2$par - z*sqrt(diag(solve(estim2$hessian))),
                    estim2$par + z*sqrt(diag(solve(estim2$hessian))))
  results <- cbind(par, rbind(conf.int[1:3, ], exp(conf.int[4:length(ini.par), ])))
  colnames(results) <- c("Estimated", "Lower", "Upper")
  
  # Generate predictions
  if (delta2) {
    beta <- par[[1]]
    gamma1 <- par[[2]]
    gamma2 <- par[[3]]
    deltau <- par[[4]]
    deltav <-  par[[5]]
    pr[vacc] <- (logistic(beta * age[vacc]) - 0.5) * (age[vacc] < kv) + 
      (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (age[vacc] - kv) + gamma2 * (age[vacc] - kv) ^ 2)) * (age[vacc] >= kv & age[vacc] < ku) + (logistic(beta * kv) - 0.5) * (0.5 + logistic(gamma1 * (ku - kv) + gamma2 * ((ku - kv) ^ 2))) * exp(- deltav * (age[vacc] - ku)) * (age[vacc] >= ku)
    pr[unvacc] <- (logistic(beta * age[unvacc]) - 0.5) * (age[unvacc] < ku) + 
      (logistic(beta * ku) - 0.5) * exp(- deltau * (age[unvacc] - ku)) * (age[unvacc] >= ku)
  } else {
    beta <- par[[1]]
    delta <- par[[2]]
    pr[vacc] <- (logistic(beta * age[vacc]) - 0.5) * (age[vacc] < kv) + 
      (logistic(beta * kv) - 0.5) * exp(- delta * (age[vacc] - kv)) * (age[vacc] >= kv)
    pr[unvacc] <- (logistic(beta * age[unvacc]) - 0.5) * (age[unvacc] < ku) + 
      (logistic(beta * ku) - 0.5) * exp(- delta * (age[unvacc] - ku)) * (age[unvacc] >= ku)
  }
  
  predicted <- pr
  return(list(results = results, convcode = estim2$convergence, lik = estim2$value, 
              y = y, ku = ku, kv = kv, 
              predicted = predicted, id = vaccinated_dummy, 
              sd = sqrt(diag(solve(estim2$hessian)))))
}

