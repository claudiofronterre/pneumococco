censorModel <- function(formula = vt_01 ~ gender + crowd_index + ses_score + year,
                       data = df, 
                       ini.par = c(0.1, 0.1, 0.5), 
                       vacc = "typ_cat2",
                       t = "t1",
                       age = "age_yrs",
                       agev = "agev",
                       tcu,
                       tcv,
                       knot = F,
                       method = "Nelder-Mead",
                       epsilon = 0.0001) {
    if(knot) {
      formula2 <- as.formula(paste0(deparse(formula, width.cutoff = 500),
                                    "+",paste(c(age, vacc), collapse= "+")))
      #excludes unvacc children aged less than tcu
      data <- data[!(data[age] < tcu & data[vacc] == 0), ] 
      data <- model.frame(formula2, data = data, na.action = na.omit)
      fit <- glm(formula, family = binomial(link = "logit"), data = data)
      D <- model.matrix(formula, data = data)
      
      p <- ncol(D)
      y <- as.numeric(model.response(data))
      id <- as.numeric(as.matrix(data[vacc]))
      age <- as.numeric(as.matrix(data[age]))
      unvacc <- which(id == 0)
      vacc <- which(id == 1)
      pr <- rep(0, nrow(data))
      logistic <- function(x) 1 / (1 + exp(- x))
      log.lik <- function(par) {
        bi <- par[1:p]
        alpha <- as.numeric(exp(D %*% bi) / (1 + exp(D %*% bi))) 
        deltau <- exp(par[p + 1])
        deltav <- exp(par[p + 2])
        beta <- par[p + 3]
        kv <- exp(par[p + 4])
        alpha2 <- exp(par[p + 5])
        
        pr[unvacc] <- alpha[unvacc] * exp(- deltau * (age[unvacc] - tcu))
        pr[vacc] <- (alpha2 + beta * age[vacc]) * (age[vacc] < kv) + 
          (alpha2 + beta * kv) * exp(- deltav * (age[vacc] - kv)) * (age[vacc] >= kv)
        
        out <- - (sum(log(1 - pr)) + sum(y * log(pr / (1 - pr))))
        as.numeric(out)
      }
      
      par.start <- c(as.numeric(fit$coefficients), log(ini.par[1]),
                     log(ini.par[2]), ini.par[3], log(ini.par[4]), log(ini.par[4]))        
      estim1 <- list(); estim1$value <- 1
      estim2 <- list(); estim2$value <- 2
      while(abs(estim2$value - estim1$value) > epsilon) {
        estim1 <- optim(par = par.start, fn = log.lik, hessian = T, control = list(maxit = 5000), 
                        method = method)  
        estim2 <- optim(par = estim1$par, fn = log.lik, hessian = T, control = list(maxit = 5000),
                        method = method)
        par.start <- estim2$par
      }
      
      
      par <- c(estim2$par[1:p], exp(estim2$par[(p + 1):(p + 2)]),
               estim2$par[p + 3], exp(estim2$par[p + 4]), exp(estim2$par[p + 5]))
      names(par) <- c(names(fit$coefficients), "deltau", "deltav", "beta", "kv", "alpha2")
      
      ## Confidence intervals 
      
      conf.level <- 0.95
      z <- - qnorm((1 - conf.level)/2)
      conf.int <- cbind(estim2$par - z*sqrt(diag(solve(estim2$hessian))),
                        estim2$par + z*sqrt(diag(solve(estim2$hessian))))
      conf.int <- rbind(conf.int[1:p, ], exp(conf.int[(p + 1):(p + 2),]),
                        conf.int[p + 3, ], 
                        exp(conf.int[(p + 4):(p + 5), ]))
      results <- cbind(par, conf.int)
      colnames(results) <- c("Estimated", "Lower", "Upper")
      rownames(results) <- c(names(fit$coefficients), "deltau", "deltav", "beta", "ku", "alpha2")
      results <- round(results, 4)
      
      bi <- as.numeric(par[1:p])
      alpha <- as.numeric(exp(D %*% bi) / (1 + exp(D %*% bi)))
      deltau <- as.numeric(par[p + 1])
      deltav <- as.numeric(par[p + 2])
      beta <-  as.numeric(par[p + 3])
      kv <- as.numeric(par[p + 4])
      alpha2 <- as.numeric(par[p + 5])
      predicted <- c(alpha[unvacc] * exp(- deltau * (age[unvacc] - tcu)),
                     (alpha2 + beta * age[vacc]) * (age[vacc] < kv) + 
                       (alpha2 + beta * kv) * exp(- deltav * (age[vacc] - kv)) * (age[vacc] >= kv))
      
      return(list(results = results, convcode = estim2$convergence, lik = estim2$value, y = y, 
                  predicted = predicted, id = id, par = estim2$par, covar = solve(estim2$hessian),
                  D = D, 
                  sd = sqrt(diag(solve(estim2$hessian)))))
    } else {
      formula2 <- as.formula(paste0(deparse(formula, width.cutoff = 500),
                                    "+",paste(c(age, vacc), collapse= "+")))
      #excludes unvacc children aged less than tcu and vacc children aged less than tcv
      data <- data[data[age] >= tcv & data[vacc] == 1 | data[age] >= tcu & data[vacc] == 0, ]
      data <- model.frame(formula2, data = data, na.action = na.omit)
      fit <- glm(formula, family = binomial(link = "logit"), data = data)
      D <- model.matrix(formula, data = data)
      
      p <- ncol(D)
      y <- as.numeric(model.response(data))
      id <- as.numeric(as.matrix(data[vacc]))
      age <- as.numeric(as.matrix(data[age]))
      unvacc <- which(id == 0)
      vacc <- which(id == 1)
      pr <- rep(0, nrow(data))
      log.lik <- function(par) {
        bi <- par[1:p]
        alpha <- as.numeric(exp(D %*% bi) / (1 + exp(D %*% bi))) 
        delta1 <- exp(par[p + 1])
        delta2 <- exp(par[p + 2])
        beta <- exp(par[p + 3]) / (1 + exp(par[p + 3]))
        
        pr[unvacc] <- alpha[unvacc] * exp(- delta1 * (age[unvacc] - tcu))
        pr[vacc] <- alpha[vacc] * beta * exp(- delta2 * (age[vacc] - tcu))
        
        out <- - (sum(log(1 - pr)) + sum(y*log(pr/(1 - pr))))
        as.numeric(out)
      }
      
      par.start <- c(as.numeric(fit$coefficients), log(ini.par[1]),
                     log(ini.par[2]), log(ini.par[3] / (1 - ini.par[3])))        
      estim1 <- list(); estim1$value <- 1
      estim2 <- list(); estim2$value <- 2
      while(abs(estim2$value - estim1$value) > epsilon) {
        estim1 <- optim(par = par.start, fn = log.lik, hessian = T, control = list(maxit = 5000), 
                        method = method)  
        estim2 <- optim(par = estim1$par, fn = log.lik, hessian = T, control = list(maxit = 5000),
                        method = method)
        par.start <- estim2$par
      }
      
      
      par <- c(estim2$par[1:p], exp(estim2$par[(p + 1):(p + 2)]),
               exp(estim2$par[p + 3])/(1 + exp(estim2$par[p + 3])))
      names(par) <- c(names(fit$coefficients), "delta1", "delta2", "beta")
      
      ## Confidence intervals 
      
      conf.level <- 0.95
      z <- - qnorm((1 - conf.level)/2)
      conf.int <- cbind(estim2$par - z*sqrt(diag(solve(estim2$hessian))),
                        estim2$par + z*sqrt(diag(solve(estim2$hessian))))
      conf.int <- rbind(conf.int[1:p, ], exp(conf.int[(p + 1):(p + 2),]),
                        exp(conf.int[p + 3, ])/(1 + exp(conf.int[p + 3, ])))
      results <- cbind(par, conf.int)
      colnames(results) <- c("Estimated", "Lower", "Upper")
      rownames(results) <- c(names(fit$coefficients), "delta1", "delta2", "beta")
      results <- round(results, 4)
      
      xi <- as.numeric(par[1:p])
      alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi)))
      delta1 <- as.numeric(par[p+1])
      delta2 <- as.numeric(par[p+2])
      beta <-  as.numeric(par[p+3])
      predicted <- c(alpha[unvacc]*exp(- delta1*(age[unvacc] - tcu)),
                     alpha[vacc]*beta*exp(- delta2*(age[vacc] - tcu)))
      
      return(list(results = results, convcode = estim2$convergence, lik = estim2$value, y = y, 
                  predicted = predicted, id = id, par = estim2$par, covar = solve(estim2$hessian),
                  D = D, 
                  sd = sqrt(diag(solve(estim2$hessian)))))
    }
}
  









