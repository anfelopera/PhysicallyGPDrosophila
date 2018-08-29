#### Utilities for physically-inspired GPs for Drosophila ####
errorMeasureRegress <- function(y, ytest, mu, varsigma,  type = "all",
                                control = list(nsigma = 1)) {
  # computing the errors
  maerror <- abs(ytest - mu)
  mserror <- (ytest - mu)^2
  std_mserror <- mserror/var(ytest)
  q2 <- sum(mserror)/sum((ytest - mean(ytest))^2)
  pva <- mserror/varsigma
  sserror <- ((ytest - mu)^2)/(2*varsigma)
  logprob <- 0.5*log(2*pi*varsigma) + sserror -
    0.5*log(2*pi*var(y)) - ((ytest - mean(y))^2)/(2*var(y))
  logprob <- logprob[complete.cases(logprob)]
  cia <- ytest >= (mu-control$nsigma*sqrt(varsigma)) &
    ytest <= (mu+control$nsigma*sqrt(varsigma))
  
  error <- c(mae = mean(maerror), mse = mean(mserror),
             smse = mean(std_mserror), msll = mean(logprob),
             q2 = 1 - q2, pva = abs(log(mean(pva))), cia = mean(cia))
  switch(type,
         mae = {return(error["mae"])},
         mse = {return(error["mse"])},
         smse = {return(error["smse"])},
         msll = {return(error["msll"])},
         q2 = {return(error["q2"])},
         pva = {return(error["pva"])},
         cia = {return(error["cia"])},
         {return(as.matrix(error, nrow = 1))}
  )
}
