
# internal function to process R2JAGS output
summarizeMPT <- function(model,
                         mcmc,
                         thetaNames,
                         covIncluded=FALSE,
                         predFactorLevels=NULL,
                         transformedParameters=NULL,
                         NgroupT1=NULL){

  summ <- mcmc$BUGSoutput$summary
  # number of parameters
  S <- max(thetaNames[,"theta"])

  idx <- ""
  if(S>1) idx <- paste0("[",1:S,"]")

  N <- dim(mcmc$BUGSoutput$mean$theta)[2]
  # which statistics to select:
  colsel <- c(1:3,5,7:9)
  uniqueNames <- thetaNames[!duplicated(thetaNames[,"theta"]),"Parameter"]

  if(covIncluded){
    selCov <- grepl("cor_", rownames(summ))
    correlation <- summ[selCov, colsel, drop=FALSE]
  }else{correlation <- NULL}

  ############################## BETA MPT SUMMARY
  if(model == "betaMPT"){
    # mean and individual parameter estimates
    mu <- summ[paste0("mean", idx),colsel, drop=FALSE]
    SD <- summ[paste0("sd",idx),colsel, drop=FALSE]
    alpha <- summ[paste0("alph",idx),colsel, drop=FALSE]
    beta <- summ[paste0("bet",idx),colsel, drop=FALSE]

    rownames(mu) <-paste0("mean_", uniqueNames)
    rownames(SD) <-paste0("sd_", uniqueNames)
    rownames(alpha) <-paste0("alpha_", uniqueNames)
    rownames(beta) <-paste0("beta_", uniqueNames)

    groupParameters <- list(mean=mu, SD=SD, alpha=alpha, beta=beta, correlation=correlation)

  }else{

    ############################## LATENT-TRAIT SUMMARY
    mu <- summ[paste0("mu",idx),colsel, drop=FALSE]
    mean <- summ[paste0("mean",idx),colsel, drop=FALSE]
    sigma <- summ[paste0("sigma",idx),colsel, drop=FALSE]
    rho <- rhoNames <- c()
    cnt <- 1
    while(cnt<S){
      rho <- rbind(rho, summ[paste0("rho[",cnt,",",(cnt+1):S,"]"),colsel, drop=FALSE])
      rhoNames <- c(rhoNames,
                    paste0("rho[",uniqueNames[cnt],",",uniqueNames[(cnt+1):S],"]"))
      cnt <- cnt+1
    }

    rownames(mu) <-paste0("latent_mu_", uniqueNames)
    rownames(mean) <-paste0("mean_", uniqueNames)
    rownames(sigma) <-paste0("latent_sigma_", uniqueNames)
    rownames(rho) <- rhoNames

    selCov <- grepl("slope_", rownames(summ))
    selFac <- grepl("factor_", rownames(summ))
    if(any(selCov)){
      slope <- summ[selCov, colsel, drop=FALSE]
    }else{slope <- NULL}
    if(any(selFac)){
      selFacSD <- grepl("SD_factor_", rownames(summ))
      factor <- summ[selFac & !selFacSD, colsel, drop=FALSE]

      # rename factor levels
      tmpNames <- sapply(rownames(factor), function(xx) strsplit(xx, split="_")[[1]][3])
      for(ff in 1:length(predFactorLevels)){
        if( !is.null(predFactorLevels[[ff]])){
          selrow <- grepl(names(predFactorLevels)[ff], tmpNames)
          rownames(factor)[selrow] <- paste0(rownames(factor)[selrow],"_", predFactorLevels[[ff]])
        }
      }


      factorSD <- summ[ selFacSD, colsel, drop=FALSE]
    }else{
      factor <- NULL ; factorSD <- NULL
    }
    rho.matrix <- mcmc$BUGSoutput$mean$rho
    dimnames(rho.matrix) <- list(uniqueNames,uniqueNames)
    groupParameters <- list(mean = mean, mu = mu, sigma = sigma, rho = rho,
                            rho.matrix=rho.matrix,
                            slope = slope, factor = factor, factorSD = factorSD, correlation=correlation)
  }
  theta.names <- apply(as.matrix(data.frame(lapply(expand.grid("theta[",1:S, ",",1:N,"]"), as.character))),
                       1, paste0, collapse="")
  individParameters <- array(data = summ[theta.names,colsel],
                             dim = c(S,N,length(colsel)))
  dimnames(individParameters) <- list(Parameter=uniqueNames,
                                      ID=1:N,
                                      Statistic=colnames(summ)[colsel])

  ############################## goodness of fit and deviance
  if(is.null(transformedParameters)){
    transPar <- NULL}
  else{
    transPar <- summ[transformedParameters,colsel, drop=FALSE]
  }

  summary <- list(groupParameters=groupParameters,
                  individParameters=individParameters,
                  fitStatistics=list(
                    "overall"=c("DIC"=mcmc$BUGSoutput$DIC,
                                "T1.observed"=mcmc$BUGSoutput$mean$T1.obs,
                                "T1.predicted"=mcmc$BUGSoutput$mean$T1.pred,
                                "p.T1"=mcmc$BUGSoutput$mean$p.T1),
                    "p.T1.individual"=mcmc$BUGSoutput$mean$p.T1ind),
                  transformedParameters=transPar)
  if(!is.null(mcmc$BUGSoutput$mean$p.T1.group)){
    summary$fitStatistics$p.T1.group <- rbind(N_per_group=NgroupT1,
                                              p.T1.group=mcmc$BUGSoutput$mean$p.T1.group)
    # colnames(summary$p.T1.individual) <- names(NgroupT1)
  }
  return(summary)
}


#' @export
print.summary.betaMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
  if(is.null(x$groupParameters)){
    warning("Clean MPT summary only available when fitting with JAGS.")
    print(x)
  }else{
    cat("Mean parameters on group level:\n")
    print(round(x$groupParameters$mean, x$round))
    cat("\nStandard deviation of parameters across individuals:\n")
    print(round(x$groupParameters$SD, x$round))
    cat("\nAlpha parameters of beta distributions:\n")
    print(round(x$groupParameters$alpha, x$round))
    cat("\nBeta parameters of beta distributions:\n\n")
    print(round(x$groupParameters$beta, x$round))

    cat("\n##############\n",
        "Overall model fit statistics (T1: Posterior predictive check):\n")
    print(round(x$fitStatistics$overall, x$round))
    cat("\nPoster predictive p-values for participants:\n")
    print(round(x$fitStatistics$p.T1.individual, x$round))
    if(!is.null(x$fitStatistics$p.T1.group)){
      cat("\nPoster predictive p-values per group:\n")
      print(round(x$fitStatistics$p.T1.group, x$round))
    }

    if(!is.null(x$transformedParameters)){
      cat("\nTransformed parameters:\n")
      print(round(x$transformedParameters, x$round))
    }
    if(!is.null(x$groupParameters$correlation) & !nrow(x$groupParameters$correlation) == 0){
      cat("\nSampled correlations of MPT parameters with covariates:\n")
      print(round(x$groupParameters$correlation, x$round))
    }

  }
}

#' @export
print.summary.traitMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
  if(is.null(x$groupParameters)){
    warning("Clean MPT summary only available when fitting with JAGS.")
    print(x)
  }else{
    cat("Mean parameters on group level:\n")
    print(round(x$groupParameters$mean, x$round))
    cat("\nMean of latent-trait values (probit-scale) across individuals:\n")
    print(round(x$groupParameters$mu, x$round))
    cat("\nStandard deviation of latent-trait values (probit scale) across individuals:\n")
    print(round(x$groupParameters$sigma, x$round))
    cat("\nCorrelations of latent-trait values on probit scale:\n")
    print(round(x$groupParameters$rho, x$round))
    cat("\nCorrelations (posterior mean estimates) in matrix form:\n")
    print(round(x$groupParameters$rho.matrix, x$round))

    cat("\n##############\n",
        "Overall model fit statistics (T1: Posterior predictive check):\n")
    print(round(x$fitStatistics$overall, x$round))
    cat("\nPoster predictive p-values for participants:\n")
    print(round(x$fitStatistics$p.T1.individual, x$round))
    if(!is.null(x$fitStatistics$p.T1.group)){
      cat("\nPoster predictive p-values per group:\n")
      print(round(x$fitStatistics$p.T1.group, x$round))
    }
    if(!is.null(x$transformedParameters)){
      cat("\nTransformed parameters:\n")
      print(round(x$transformedParameters, x$round))
    }

    if(!is.null(x$groupParameters$slope)){
      cat("\nSlope parameters for predictor variables:\n")
      print(round(x$groupParameters$slope, x$round))
    }

    if(!is.null(x$groupParameters$factor)){
      cat("\nEffects of factors on latent scale (additive shift from overall mean):\n")
      print(round(x$groupParameters$factor, x$round))
      cat("\nFactor SD on latent scale:\n")
      print(round(x$groupParameters$factorSD, x$round))
    }

    if(!is.null(x$groupParameters$correlation) & !nrow(x$groupParameters$correlation) == 0){
      cat("\nSampled correlations of covariates with MPT parameters:\n")
      print(round(x$groupParameters$correlation, x$round))
    }
  }
}

#' @export
summary.betaMPT <- function(object, round=3, ...){
#   if(!object$sampler %in% c("jags", "JAGS")){
#     warning("Clean MPT summary only available when fitting with JAGS.")
#   }
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  class(summ) <- "summary.betaMPT"
  return(summ)
}

#' @export
summary.traitMPT <- function(object, round=3, ...){
#   if(!object$sampler %in% c("jags", "JAGS")){
#     warning("Clean MPT summary only available when fitting with JAGS.")
#   }
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  class(summ) <- "summary.traitMPT"
  return(summ)
}

#' @export
print.betaMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
#   if(!x$sampler %in% c("jags", "JAGS")){
#     warning("\nClean MPT summary only available when fitting with JAGS.")
#   }else{
    print(round(cbind("Group Mean" = x$summary$groupParameters$mean[,1],
                      "Group SD" = x$summary$groupParameters$SD[,1]),4))

  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}

#' @export
print.traitMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
#   if(!x$sampler %in% c("jags", "JAGS")){
#     warning("\nClean MPT summary only available when fitting with JAGS.")
#   }else{
    print(round(cbind("Mean(MPT Parameters)" = x$summary$groupParameters$mean[,1],
                      "Mu(latent-traits)" = x$summary$groupParameters$mu[,1],
                      "SD(latent-traits)" = x$summary$groupParameters$sigma[,1]),4))

  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}
