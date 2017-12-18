#' Summarize JAGS Output for Hierarchical MPT Models
#'
#' Provide clean and readable summary statistics tailored to MPT models based on the JAGS output.
#'
# @param model either \code{"betaMPT"} or \code{"traitMPT"}
#' @inheritParams summarizeMCMC
#' @param mcmc the actual mcmc.list output of the sampler of a fitted MPT model (accesible via \code{fittedModel$runjags$mcmc})
#' @param mptInfo the internally stored information about the fitted MPT model (accesible via \code{fittedModel$mptInfo})
# @param dic whether to compute DIC statistic for model selection (requires additional sampling!)
# @param M number of posterior predictive samples to compute T1 statistic
#' @param summ optional argument for internal use
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
#' @import rjags
summarizeMPT <- function(mcmc, mptInfo, probs = c(.025,.50,.975), summ = NULL){
  if(is.null(summ) | !all(paste0(probs*100, "%") %in% colnames(summ)))
    summ <- summarizeMCMC(mcmc, probs = probs)

  predFactorLevels <- mptInfo$predFactorLevels
  transformedParameters <- mptInfo$transformedParameters
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
  uniqueNames <- mptInfo$thetaUnique

  selCov <- grep("cor_", rownames(summ))
  if(length(selCov) != 0){
    correlation <- summ[selCov, , drop=FALSE]
  }else{correlation <- NULL}

  selFE <- grep("thetaFE", rownames(summ))
  if(length(selFE) != 0){
    thetaFE <- summ[selFE, , drop=FALSE]
    if(length(mptInfo$thetaFixed) == nrow(thetaFE)){
      rownames(thetaFE) <- paste0("thetaFE_", mptInfo$thetaFixed)
    }
  }else{thetaFE <- NULL}

  rho <- rhoNames <- c()
  cnt <- 1
  while(cnt<S){
    rho <- rbind(rho, summ[paste0("rho[",cnt,",",(cnt+1):S,"]"),, drop=FALSE])
    rhoNames <- c(rhoNames,
                  paste0("rho[",uniqueNames[cnt],",",uniqueNames[(cnt+1):S],"]"))
    cnt <- cnt+1
  }
  if(S == 1) rho <- matrix(1,1,1, dimnames=list(uniqueNames,uniqueNames))
  rownames(rho) <- rhoNames
  rho.matrix <- getRhoMatrix(uniqueNames, rho)

  mean <- summ[paste0("mean", idx),, drop=FALSE]
  rownames(mean) <-paste0("mean_", uniqueNames)
  ############################## BETA MPT SUMMARY
  if(model == "simpleMPT"){
    SD <- summ[paste0("sd",idx),, drop=FALSE]
    rownames(SD) <-paste0("sd_", uniqueNames)
    groupParameters <- list(mean=mean, SD=SD,
                            rho=rho, rho.matrix = rho.matrix,
                            correlation=correlation, thetaFE=NULL)

  }else if(model == "betaMPT"){
    # mean and individual parameter estimates
    SD <- summ[paste0("sd",idx),, drop=FALSE]
    alpha <- summ[paste0("alph",idx),, drop=FALSE]
    beta <- summ[paste0("bet",idx),, drop=FALSE]

    rownames(SD) <-paste0("sd_", uniqueNames)
    rownames(alpha) <-paste0("alpha_", uniqueNames)
    rownames(beta) <-paste0("beta_", uniqueNames)

    # rho <- summgetRhoBeta(uniqueNames, mcmc, corProbit=mptInfo$corProbit)
    # rho.matrix <- getRhoMatrix(uniqueNames, rho)

    groupParameters <- list(mean=mean, SD=SD, alpha=alpha, beta=beta,
                            rho=rho, rho.matrix = rho.matrix,
                            correlation=correlation, thetaFE=thetaFE)

  }else if (model == "traitMPT"){

    ############################## LATENT-TRAIT SUMMARY
    mu <- summ[paste0("mu",idx),, drop=FALSE]
    sigma <- summ[paste0("sigma",idx),, drop=FALSE]

    rownames(mu) <-paste0("latent_mu_", uniqueNames)
    rownames(sigma) <-paste0("latent_sigma_", uniqueNames)

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
          selrow <- grepl(paste0(names(predFactorLevels)[ff],"["), tmpNames, fixed=TRUE)
          rownames(factor)[selrow] <- paste0(rownames(factor)[selrow],"_", predFactorLevels[[ff]])
        }
      }


      factorSD <- summ[ selFacSD, , drop=FALSE]
    }else{
      factor <- NULL ; factorSD <- NULL
    }

    groupParameters <- list(mean = mean, mu = mu, sigma = sigma, rho = rho,
                            rho.matrix=rho.matrix,
                            slope = slope, factor = factor, factorSD = factorSD,
                            correlation=correlation, thetaFE=thetaFE)
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
                  fitStatistics=NULL,
                  transformedParameters=transPar)

  summary$call <- "(summarizeMPT called manually)"
  summary$round <- 3
  class(summary) <- paste0("summary.", model)

  return(summary)
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

