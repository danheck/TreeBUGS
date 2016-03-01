#' Summarize JAGS Output for Hierarchical MPT Models
#'
#' Provide clean and readable summary statistics tailored to MPT models based on the JAGS output.
#'
# @param model either \code{"betaMPT"} or \code{"traitMPT"}
#' @param mcmc the actual mcmc.list output of the sampler of a fitted MPT model (accesible via \code{fittedModel$runjags$mcmc})
#' @param mptInfo the internally stored information about the fitted MPT model (accesible via \code{fittedModel$mptInfo})
# @param dic whether to compute DIC statistic for model selection (requires additional sampling!)
#' @param summ optional for internal use
# @param ... further arguments passed to \code{\link[coda]{dic.samples}} (e.g., \code{n.iter})
#' @details The MPT-specific summary is computed directly after fitting a model. However, this function might be used manually after removing MCMC samples (e.g., extending the burnin period).
#' @examples
#' # Remove additional burnin samples and recompute MPT summary
#' \dontrun{
#' # start later or thin (see ?window)
#' mcmc.subsamp <- window(fittedModel$runjags$mcmc, start=3001, thin=2)
#' new.mpt.summary <- summarizeMPT(mcmc.subsamp, fittedModel$mptInfo)
#' new.mpt.summary
#' }
#' @export
summarizeMPT <- function(mcmc,
                         mptInfo,
                         summ = NULL
                         ){
  if(is.null(summ))
    summ <- summarizeMCMC(mcmc)

  predFactorLevels <- mptInfo$predFactorLevels
  transformedParameters <- mptInfo$transformedParameters
  NgroupT1 <- mptInfo$T1group$NgroupT1
  thetaNames <- mptInfo$thetaNames
  model <- mptInfo$model

  # summary <- summary(mcmc) # mcmc$BUGSoutput$summary
  # number of parameters
  S <- max(thetaNames[,"theta"])

  idx <- ""
  if(S>1) idx <- paste0("[",1:S,"]")

  N <- nrow(mptInfo$data)
  # which statistics to select:
  # colsel <- c(1:3,5,7:9) # for R2JAGS
  # colsel <- c("Mean","SD","Lower95","Median","Upper95","SSeff","psrf") # for runjags
  uniqueNames <- mptInfo$thetaUnique # thetaNames[!duplicated(thetaNames[,"theta"]),"Parameter"]

  selCov <- grep("cor_", rownames(summ))
  if(length(selCov) != 0){
    correlation <- summ[selCov, , drop=FALSE]
  }else{correlation <- NULL}

  ############################## BETA MPT SUMMARY
  if(model == "betaMPT"){
    # mean and individual parameter estimates
    mean <- summ[paste0("mean", idx),, drop=FALSE]
    SD <- summ[paste0("sd",idx),, drop=FALSE]
    alpha <- summ[paste0("alph",idx),, drop=FALSE]
    beta <- summ[paste0("bet",idx),, drop=FALSE]

    rownames(mean) <-paste0("mean_", uniqueNames)
    rownames(SD) <-paste0("sd_", uniqueNames)
    rownames(alpha) <-paste0("alpha_", uniqueNames)
    rownames(beta) <-paste0("beta_", uniqueNames)

    groupParameters <- list(mean=mean, SD=SD, alpha=alpha, beta=beta, correlation=correlation)

  }else{

    ############################## LATENT-TRAIT SUMMARY
    mu <- summ[paste0("mu",idx),, drop=FALSE]
    mean <- summ[paste0("mean",idx),, drop=FALSE]
    sigma <- summ[paste0("sigma",idx),, drop=FALSE]
    rho <- rhoNames <- c()
    cnt <- 1
    while(cnt<S){
      rho <- rbind(rho, summ[paste0("rho[",cnt,",",(cnt+1):S,"]"),, drop=FALSE])
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
      slope <- summ[selCov, , drop=FALSE]
    }else{slope <- NULL}
    if(any(selFac)){
      selFacSD <- grepl("SD_factor_", rownames(summ))
      factor <- summ[selFac & !selFacSD, , drop=FALSE]

      # rename factor levels
      tmpNames <- sapply(rownames(factor), function(xx) strsplit(xx, split="_")[[1]][3])
      for(ff in 1:length(predFactorLevels)){
        if( !is.null(predFactorLevels[[ff]])){
          selrow <- grepl(names(predFactorLevels)[ff], tmpNames)
          rownames(factor)[selrow] <- paste0(rownames(factor)[selrow],"_", predFactorLevels[[ff]])
        }
      }


      factorSD <- summ[ selFacSD, , drop=FALSE]
    }else{
      factor <- NULL ; factorSD <- NULL
    }

    rho.matrix <- getRhoMatrix(uniqueNames, rho)
    groupParameters <- list(mean = mean, mu = mu, sigma = sigma, rho = rho,
                            rho.matrix=rho.matrix,
                            slope = slope, factor = factor, factorSD = factorSD, correlation=correlation)
  }
  theta.names <- apply(as.matrix(data.frame(lapply(expand.grid("theta[",1:S, ",",1:N,"]"), as.character))),
                       1, paste0, collapse="")
  individParameters <- array(data = summ[theta.names,],
                             dim = c(S,N,ncol(summ)))
  dimnames(individParameters) <- list(Parameter=uniqueNames,
                                      ID=1:N,
                                      Statistic=colnames(summ))

  ############################## goodness of fit and deviance
  if(is.null(transformedParameters)){
    transPar <- NULL
  }else{
    transPar <- summ[transformedParameters,, drop=FALSE]
  }
  # dic <- extract(mcmc, "dic")
  summary <- list(groupParameters=groupParameters,
                  individParameters=individParameters,
                  fitStatistics=list(
                    "overall"=c(
                      # "DIC"=sum(dic$deviance) + mean( sum(dic$penalty)), #$BUGSoutput$DIC,
                      "T1.observed"=summ["T1.obs","Mean"],
                      "T1.predicted"=summ["T1.pred","Mean"],
                      "p.T1"=summ["p.T1","Mean"]),
                    "p.T1.individual"=summ[paste0("p.T1ind[",1:N,"]"),"Mean"], #mcmc$BUGSoutput$mean$p.T1ind),
                  transformedParameters=transPar))
  selT1group <- grep("p.T1.group", rownames(summ))
  if(length(selT1group) != 0){
  # if(!is.null(mcmc$BUGSoutput$mean$p.T1.group)){
    summary$fitStatistics$p.T1.group <- rbind(N_per_group=NgroupT1,
                                              p.T1.group=summ[selT1group,"Mean"])
    # colnames(summary$p.T1.individual) <- names(NgroupT1)
  }

  summary$call <- "(summarizeMPT called manually)"
  summary$round <- 3
  class(summary) <- paste0("summary.", model)

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
    if(!is.null(x$dic)){
      print(x$dic)
    }
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
    if(!is.null(x$groupParameters$correlation) && !nrow(x$groupParameters$correlation) == 0){
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
    if(!is.null(x$dic)){
      print(x$dic)
    }
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

    if(!is.null(x$groupParameters$correlation) && !nrow(x$groupParameters$correlation) == 0){
      cat("\nSampled correlations of covariates with MPT parameters:\n")
      print(round(x$groupParameters$correlation, x$round))
    }
  }
}

#' @export
summary.betaMPT <- function(object, round=3, ...){
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  # class(summ) <- "summary.betaMPT"
  return(summ)
}

#' @export
summary.traitMPT <- function(object, round=3, ...){
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
#   class(summ) <- "summary.traitMPT"
  return(summ)
}

#' @export
print.betaMPT <- function(x,  ...){
  cat("Call: \n")
  print(x$call)
  cat("\n")
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





getRhoMatrix <- function (uniqueNames, rho) {
  S <- length(uniqueNames)
  rho.matrix <- matrix(1, S, S, dimnames=list(uniqueNames,uniqueNames))
  if(S>1){
    cnt <- 0
    for(s1 in 1:(S-1)){
      for(s2 in (s1+1):S){
        rho.matrix[s1,s2] <- rho.matrix[s2,s1]  <- rho[cnt <- cnt+1]
      }
    }
  }
  rho.matrix
}


# own MCMC summary
summarizeMCMC <- function(mcmc){
  # summ <- summary(mcmc, quantiles = c(0.025, 0.5, 0.975))
  mcmc.mat <- do.call("rbind", mcmc)
  summTab <- cbind("Mean"=apply(mcmc.mat,2,mean),
                   "SD"=apply(mcmc.mat,2,sd),
                   t(apply(mcmc.mat, 2, quantile, c(.05,.5,.95))),
                   "Time-series SE"=NA, "n.eff" = NA ,
                   "Rhat" = NA, "R_95%"=NA)
  #     summ[[1]][,1:2], summ[[2]], "Time-series SE"=summ[[1]][,4]
#   rm(summ)
#   gc(verbose=FALSE)
  rn <- rownames(summTab)
  sel.notT1 <- setdiff(1:nrow(summTab), union(grep("T1", rn), grep(".pred.mean", rn)))
  try({
    summTab[sel.notT1,7] <- round(effectiveSize(mcmc[,sel.notT1]))
    summTab[sel.notT1,6] <- summTab[sel.notT1,2] / sqrt(summTab[sel.notT1,7]  )
  })
  try( summTab[sel.notT1,8:9] <- gelman.diag(mcmc[,sel.notT1], multivariate=FALSE)[[1]])
#   if(any(is.na(summTab[,"Rhat"])))
#     warning("Gelman-Rubin convergence diagnostic Rhat could not be computed.")
  # n.eff <-

  summTab
}
