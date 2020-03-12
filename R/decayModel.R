decayModel <- function(formula = vt_01 ~ gender + crowd_index + ses_score + year,
                       data = df, 
                       ini.par = c(0.1, 0.1, 0.5), 
                       vacc = "typ_cat2",
                       t = "t1",
                       age = "age_yrs",
                       agev = "agev",
                       tc = 0.5,
                       gamma,
                       delta2,
                       method = "Nelder-Mead",
                       epsilon = 0.0001,
                       eig = F) {
  
if(gamma == FALSE & delta2 == FALSE) {
  formula2 <- as.formula(paste0(deparse(formula, width.cutoff = 500),
                                "+",paste(c(age, vacc), collapse= "+")))
  data <- model.frame(formula2, data = data)
  data <- data[data[age] >= tc, ] #excludes children aged < tc
  fit <- glm(formula, family = binomial(link = "logit"), data = data)
  D <- model.matrix(formula, data = data)
  
  if(eig == T) {
    ev <- eigen(var(D))$values
    plot(cumsum(ev)/sum(ev))
    print(kappa(var(D)))
  } 
  
  p <- ncol(D)
  y <- as.numeric(model.response(data))
  id <- as.numeric(as.matrix(data[vacc]))
  age <- as.numeric(as.matrix(data[age]))
  
log.lik <- function(par) {
   xi <- par[1:p]
   alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi)))
   delta <- exp(par[p+1])
   beta <- exp(par[p+2])/(1 + exp(par[p+2]))

   pr <- alpha*(beta^id)*exp(-delta*(age - tc))

   out <- - (sum(log(1 - pr)) + sum(y*log(pr/(1 - pr))))
   as.numeric(out)
  }

par.start <- c(as.numeric(fit$coefficients), 
               log(ini.par[1]),
               log(ini.par[2]/(1 - ini.par[2])))        

estim1 <- list(); estim1$value <- 1
estim2 <- list(); estim2$value <- 2

while(abs(estim2$value - estim1$value) > epsilon) {
    estim1 <- optim(par = par.start, fn = log.lik, hessian = T, control = list(maxit = 5000), method = method)  
    estim2 <- optim(par = estim1$par, fn = log.lik, hessian = T, control = list(maxit = 5000), method = method)
    par.start <- estim2$par
  }
  
  
  par <- c(estim2$par[1:p], exp(estim2$par[p + 1]),
           exp(estim2$par[p + 2])/(1 + exp(estim2$par[p + 2])))
  names(par) <- c(names(fit$coefficients), "delta", "beta")
  
  ## Confidence intervals 
  
  conf.level <- 0.95
  z <- - qnorm((1 - conf.level)/2)
  conf.int <- cbind(estim2$par - z*sqrt(diag(solve(estim2$hessian))),
                    estim2$par + z*sqrt(diag(solve(estim2$hessian))))
  conf.int <- rbind(conf.int[1:p, ], exp(conf.int[(p + 1),]),
                    exp(conf.int[p + 2, ])/(1 + exp(conf.int[p + 2, ])))
  results <- cbind(par, conf.int)
  colnames(results) <- c("Estimated", "Lower", "Upper")
  rownames(results) <- c(names(fit$coefficients), "delta", "beta")
  results <- round(results, 4)
  
  xi <- as.numeric(par[1:p])
  alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi)))
  delta <- as.numeric(par[p+1])
  beta <-  as.numeric(par[p+2])
  predicted <- alpha*beta^id*exp(-delta*age)
  
  return(list(results = results, convcode = estim2$convergence, lik = estim2$value, y = y, 
              predicted = predicted, id = id, covar = solve(estim2$hessian),
              sd = sqrt(diag(solve(estim2$hessian)))))
  
  
} else {
  if (delta2 == FALSE & gamma == T) {
    if (sum(is.na(data[t])  & data[vacc] == 1) > 0) {
      data <- data[-which(is.na(data[t]) & data[vacc] == 1), ] #remove vacc children without age at vaccination
    }
    data[t][is.na(data[t])] <- 0 
    data[agev][is.na(data[agev])] <- 0
    formula2 <- as.formula(paste0(deparse(formula, width.cutoff = 500),
                                  "+",paste(c(t, age, vacc, agev), collapse= "+")))
    data <- model.frame(formula2, data = data)
    data <- data[data[age] >= tc & data[agev] <= tc & data[agev] >= 0, ] #excludes children aged < tc and vacc afetr tc
    fit <- glm(formula, family = binomial(link = "logit"), data = data[data[vacc] == 0, ])
    D <- model.matrix(formula, data = data)
    
    if(eig == T) {
      ev <- eigen(var(D))$values
      plot(cumsum(ev)/sum(ev))
      print(kappa(var(D)))
    } 
    
    p <- ncol(D)
    y <- as.numeric(model.response(data)) - 1
    id <- as.numeric(as.matrix(data[vacc]))
    age <- as.numeric(as.matrix(data[age]))
    t <- as.numeric(as.matrix(data[t]/365))
    
    log.lik <- function(par) {
      xi <- par[1:p]
      alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi))) 
      delta <- exp(par[p+1])
      gamma <- exp(par[p+2])
      beta <- exp(par[p+3])/(1 + exp(par[p+3]))
      
      p <- alpha*(beta^id)*exp(-delta*(age - tc) -gamma*t*id)
      
      out <- - (sum(log(1 - p)) + sum(y*log(p/(1 - p))))
      as.numeric(out)
    }
    
    par.start <- c(as.numeric(fit$coefficients), log(ini.par[1]),
                   log(ini.par[2]), log(ini.par[3]/(1 - ini.par[3])))        
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
    names(par) <- c(names(fit$coefficients), "delta", "gamma", "beta")
    
    ## Confidence intervals 
    
    conf.level <- 0.95
    z <- - qnorm((1 - conf.level)/2)
    conf.int <- cbind(estim2$par - z*sqrt(diag(solve(estim2$hessian))),
                      estim2$par + z*sqrt(diag(solve(estim2$hessian))))
    conf.int <- rbind(conf.int[1:p, ], exp(conf.int[(p + 1):(p + 2),]),
                      exp(conf.int[p + 3, ])/(1 + exp(conf.int[p + 3, ])))
    results <- cbind(par, conf.int)
    colnames(results) <- c("Estimated", "Lower", "Upper")
    rownames(results) <- c(names(fit$coefficients), "delta", "gamma", "beta")
    results <- round(results, 4)
    
    xi <- as.numeric(par[1:p])
    alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi)))
    delta <- as.numeric(par[p+1])
    beta <-  as.numeric(par[p+2])
    predicted <- alpha*beta^id*exp(-delta*age -gamma*t*id)
    
    return(list(results = results, convcode = estim2$convergence, lik = estim2$value, y = y, 
                predicted = predicted, id = id, sd = sqrt(diag(solve(estim2$hessian)))))
    
  } else {
    formula2 <- as.formula(paste0(deparse(formula, width.cutoff = 500),
                                  "+",paste(c(age, vacc), collapse= "+")))
    data <- data[data[age] >= tc, ] #excludes children aged < tc
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
      xi <- par[1:p]
      alpha <- as.numeric(exp(D%*%xi)/(1 + exp(D%*%xi))) 
      delta1 <- exp(par[p+1])
      delta2 <- exp(par[p+2])
      beta <- exp(par[p+3])/(1 + exp(par[p+3]))
      
      pr[unvacc] <- alpha[unvacc] * exp(- delta1 * (age[unvacc] - tc))
      pr[vacc] <- alpha[vacc] * beta * exp(- delta2 * (age[vacc] - tc))
      
      out <- - (sum(log(1 - pr)) + sum(y*log(pr/(1 - pr))))
      as.numeric(out)
    }
    
    par.start <- c(as.numeric(fit$coefficients), log(ini.par[1]),
                   log(ini.par[2]), log(ini.par[3]/(1 - ini.par[3])))        
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
    predicted <- c(alpha[unvacc]*exp(- delta1*(age[unvacc] - tc)),
                   alpha[vacc]*beta*exp(- delta2*(age[vacc] - tc)))
    
    return(list(results = results, convcode = estim2$convergence, lik = estim2$value, y = y, 
                predicted = predicted, id = id, par = estim2$par, covar = solve(estim2$hessian),
                D = D, 
                sd = sqrt(diag(solve(estim2$hessian)))))
    
  }
  
}
  }









