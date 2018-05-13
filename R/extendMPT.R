#' Extend MCMC Sampling for MPT Model
#'
#' Adds more MCMC samples to the fitted MPT model.
#'
#' @param fittedModel a fitted \code{\link{traitMPT}} or \code{\link{betaMPT}}
#' @inheritParams betaMPT
#' @param ... further arguments passed to \code{extend.jags} (see arguments listed in: \link[runjags]{run.jags}).
#'
#' When drawing more samples, JAGS requires an additional adaptation phase, in which the MCMC
#' sampling procedure is adjusted. Note that the MCMC sampling will still
#' give correct results even if the warning appears: "Adaptation incomplete."
#' (this just means that sampling efficiency is not optimal).
#'
#' @export
extendMPT <- function(fittedModel, n.iter = 10000, n.adapt = 1000, n.burnin = 0, ...){

  args <- list(...)
  if ("n.thin" %in% names(args))
    warning("Thinnning interval cannot be changed and is ignored!")
  args$n.thin <- args$burnin <- args$adapt <- args$sample <- args$summarise <- args$model <- NULL

  # remove correlations (otherwise, extension not possible)
  sel.cor <- grep("cor_", varnames(fittedModel$runjags$mcmc), fixed=TRUE)
  if(class(fittedModel) == "betaMPT")
    sel.cor <- c(sel.cor, grep("rho", varnames(fittedModel$runjags$mcmc), fixed=TRUE))
  if(length(sel.cor)>0)
    fittedModel$runjags$mcmc <- fittedModel$runjags$mcmc[,- sel.cor]

  args_extend <- c(list(runjags.object = fittedModel$runjags,
                        burnin = n.burnin,
                        sample = ceiling((n.iter-n.burnin)/fit$runjags$thin),
                        adapt = n.adapt,
                        summarise = FALSE),
                   args)
  fittedModel$runjags <- do.call("extend.jags", args_extend)

  # add correlations
  covData <- fittedModel$mptInfo$covData
  predTable <- fittedModel$mptInfo$predTable
  if(!is.null(covData) | fittedModel$mptInfo$model == "betaMPT"){
    if(!is.null(predTable) & fittedModel$mptInfo$model == "traitMPT"){
      isPred <- (1:ncol(covData)) %in% predTable$covIdx
    }else{
      isPred <- rep(FALSE, length(fittedModel$mptInfo$predType))
    }

    sel <- fittedModel$mptInfo$predType == "c" & !isPred
    if(any(sel) | class(fittedModel) == "betaMPT"){
      cdat <- covData[,sel,drop = FALSE]
      fittedModel$runjags$mcmc <- as.mcmc.list(
        lapply(fittedModel$runjags$mcmc, corSamples,
               covData=cdat,
               thetaUnique=fittedModel$mptInfo$thetaUnique,
               rho=ifelse(class(fittedModel) == "betaMPT", TRUE, FALSE),
               corProbit = fittedModel$mptInfo$corProbit))
    }
  }

  fittedModel$mcmc.summ <- summarizeMCMC(fittedModel$runjags$mcmc)
  fittedModel$summary <- summarizeMPT(mcmc = fittedModel$runjags$mcmc,
                                      summ = fittedModel$mcmc.summ,
                                      mptInfo = fittedModel$mptInfo)
  fittedModel$call <- c(fittedModel$call, match.call())
  fittedModel
}

