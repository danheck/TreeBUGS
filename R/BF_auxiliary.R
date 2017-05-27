
# Approximate posterior distribution of MPT parameters by a well-known, simple density
# Curently, only distribution="beta" is implemented
approximatePosterior <- function(mod,
                                 sample=2000,
                                 distribution="beta",
                                 lower=.1,
                                 upper=1000){
  # estimate alpha/beta parameters of beta approximation
  S <- length(mod$mptInfo$thetaUnique)
  betapar <- matrix(1, S, 2, dimnames=list(mod$mptInfo$thetaUnique,
                                           c("alpha","beta")))
  for(i in 1:S){
    sel <- paste0("theta[",i,",1]")
    ss <- unlist(mod$runjags$mcmc[,sel])
    ss <- sample(ss, min(sample, length(ss)))
    suppressWarnings(betapar[i,] <- fitdistr(ss,"beta",
                                             list(shape1=1,shape2=1),
                                             lower=lower,
                                             upper=upper)$estimate)
  }
  betapar
}


# resample MCMC iterations from simpleMPT object
resampling <- function(mod,
                       resample=1000){
  S <- length(mod$mptInfo$thetaUnique)

  sel <- paste0("theta[",1:S, ",1]")
  C <- length(mod$runjags$mcmc)
  R <- nrow(mod$runjags$mcmc[[1]])
  r <- ceiling(resample/C)
  if(resample > C*R){
    warning("Fitted models have less samples than required for resampling.",
            "Posterior samples will be reused!")
    rr <- lapply(mod$runjags$mcmc[,sel, drop=FALSE],
                 function(mm) mm[sample(1:R, r, replace=TRUE),,drop=FALSE])
  }else{
    rr <- lapply(mod$runjags$mcmc[,sel, drop=FALSE],
                 function(mm) mm[sample(1:R, r),,drop=FALSE])
  }
  do.call("rbind", rr)[1:resample,,drop=FALSE]
}

