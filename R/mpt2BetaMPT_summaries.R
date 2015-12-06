
# internal function to process R2JAGS output
summarizeBetaMPT <- function(mcmc, thetaNames, sampler="JAGS", transformedParameters=NULL){

  # number of parameters
  S <- max(thetaNames[,"theta"])

  if(sampler %in% c("JAGS","jags")){
    N <- dim(mcmc$BUGSoutput$mean$theta)[2]
    # which statistics to select:
    colsel <- c(1:3,5,7:9)

    # mean and individual parameter estimates
    mu <- mcmc$BUGSoutput$summary[paste0("mnb[",1:S,"]"),colsel]
    variance <- mcmc$BUGSoutput$summary[paste0("varp[",1:S,"]"),colsel]
    alpha <- mcmc$BUGSoutput$summary[paste0("alph[",1:S,"]"),colsel]
    beta <- mcmc$BUGSoutput$summary[paste0("bet[",1:S,"]"),colsel]

    uniqueNames <- thetaNames[!duplicated(thetaNames[,"theta"]),"Par"]
    rownames(mu) <-paste0("mean_", uniqueNames)
    rownames(variance) <-paste0("var_", uniqueNames)
    rownames(alpha) <-paste0("alpha_", uniqueNames)
    rownames(beta) <-paste0("beta_", uniqueNames)
    meanParameters <- list(mean=mu, variance=variance, alpha=alpha, beta=beta)

    individParameters <- array(c(mcmc$BUGSoutput$mean$theta,
                                 mcmc$BUGSoutput$median$theta,
                                 mcmc$BUGSoutput$sd$theta), c(S,N,3))
    dimnames(individParameters) <- list(Parameter=uniqueNames,
                                        ID=1:N, Statistic=c("Mean","Median","SD"))

    # goodness of fit and deviance
    if(is.null(transformedParameters)){
      transPar <- NULL}
    else{
      transPar <- mcmc$BUGSoutput$summary[transformedParameters,colsel, drop=FALSE]
    }

    summary <- list(meanParameters=meanParameters,
                    individParameters=individParameters,
                    fitStatistics=list(
                      "overall"=c("DIC"=mcmc$BUGSoutput$DIC,
                                  "T1.observed"=mcmc$BUGSoutput$mean$T1.obs,
                                  "T1.predicted"=mcmc$BUGSoutput$mean$T1.pred,
                                  "p.T1"=mcmc$BUGSoutput$mean$p.T1),
                      "p.T1individual"=mcmc$BUGSoutput$mean$p.T1ind),
                    transformedParameters=transPar)
  }else{
    warning("Clean MPT summary statistics only available when using JAGS.")
    summary <- mcmc$BUGSoutput$summary
  }
  return(summary)
}


#' @export
print.summary.betaMPT <- function(x,  ...){
  if(is.null(x$meanParameters)){
    warning("Clean MPT summary only available when fitting with JAGS.")
    print(x)
  }else{
    cat("Mean parameters on group level:\n")
    print(round(x$meanParameters$mean, 5))
    cat("\nVariance of parameters across individuals:\n")
    print(round(x$meanParameters$variance, 5))
    cat("\nOverall model fit statistics (T1: Posterior predictive check):\n")
    print(round(x$fitStatistics$overall, 4))
    cat("\nPoster predictive p-values for individual data sets:\n")
    print(round(x$fitStatistics$p.T1individual, 4))
    if(!is.null(x$transformedParameters)){
      cat("\nTransformed parameters:\n")
      print(round(x$transformedParameters, 4))
    }

  }
}

#' @export
summary.betaMPT <- function(object,  ...){
  if(!object$sampler %in% c("jags", "JAGS")){
    warning("Clean MPT summary only available when fitting with JAGS.")
  }
  summ <- object$summary
  class(summ) <- "summary.betaMPT"
  return(summ)
}

#' @export
print.betaMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
  if(!x$sampler %in% c("jags", "JAGS")){
    warning("\nClean MPT summary only available when fitting with JAGS.")
  }else{
    print(round(cbind("Group Mean" = x$summary$meanParameters$mean[,1],
                "Group Variance" = x$summary$meanParameters$variance[,1]),5))
  }
  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}
