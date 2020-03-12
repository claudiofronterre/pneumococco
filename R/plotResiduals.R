plotResiduals <- function(model, ngroups, size, legend) {
  
  y <- model$y
  phat <- model$predicted
  id <- as.logical(model$id)
  
  vacc <- data.frame(phat = phat[id], y = y[id])
  unvacc <- data.frame(phat = phat[!id], y = y[!id])
  vacc <- vacc[order(vacc[,1]), ]  
  unvacc <- unvacc[order(unvacc[,1]), ]  
  
  library(zoo)
  
  #Vaccinated
  k <- floor(nrow(vacc) / ngroups)
  ysum <- rollapply(vacc$y, width = k,  FUN = sum, by = k)
  psum.vacc <- rollapply(vacc$phat, width = k, FUN = sum, by = k)
  den <- rollapply(vacc$phat * (1 - vacc$phat), width = k, FUN = sum, by = k)
  resid.vacc <- (ysum - psum.vacc) / sqrt(den)
  
  #Unvaccinated
  k <- floor(nrow(unvacc) / ngroups)
  ysum <- rollapply(unvacc$y, width = k,  FUN = sum, by = k)
  psum.unvacc <- rollapply(unvacc$phat, width = k, FUN = sum, by = k)
  den <- rollapply(unvacc$phat*(1 - unvacc$phat), width = k, FUN = sum, by = k)
  resid.unvacc <- (ysum - psum.unvacc) / sqrt(den)
  
  res <- data.frame(residuals = c(resid.vacc,resid.unvacc), 
                    predicted = c(psum.vacc, psum.unvacc),
                    group =  c(rep("vaccinated", length(resid.vacc)), 
                               rep("unvaccinated", length(resid.unvacc))))
  
  
  library(ggpubr)
  pal <- get_palette(palette = "lancet", 2)
  m <- max(abs(res$residuals))
  yl <- ifelse(m < floor(m) + 0.5, ceiling(m) - 0.5, ceiling(m))
  ggplot(data = res, aes(x = predicted, y = residuals, color = group)) + 
    geom_point(size = 1.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    ylim(c(-yl, yl)) +
    scale_color_manual(values = rev(pal), name = "", lab = c("Unvaccinated", "Vaccinated")) +
    theme_bw() +
    theme(legend.position = legend) 
    
}

