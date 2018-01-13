
#' Extend MCMC Sampling for MPT Model
#'
#' @param fittedModel a fitted \code{\link{traitMPT}} or \code{\link{betaMPT}}
#' @inheritParams betaMPT
#' @param ... further arguments passed to \link[runjags]{run.jags}
#' @export
extendMPT <- function(fittedModel, n.iter = 10000, n.adapt = 1000,
                      n.burnin = 0, n.thin = 5, ...){

  # remove correlations (otherwise, extension not possible)
  sel.cor <- grep("cor_", varnames(fittedModel$runjags$mcmc), fixed=TRUE)
  if(class(fittedModel) == "betaMPT")
    sel.cor <- c(sel.cor, grep("rho", varnames(fittedModel$runjags$mcmc), fixed=TRUE))
  if(length(sel.cor)>0)
    fittedModel$runjags$mcmc <- fittedModel$runjags$mcmc[,- sel.cor]
  tmp <- extend.jags(fittedModel$runjags,
                     burnin = n.burnin,
                     sample = ceiling((n.iter-n.burnin)/n.thin),
                     adapt = n.adapt,
                     thin=n.thin,
                     summarise = FALSE, ...)
  fittedModel$runjags <- tmp

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
      tmp <- as.mcmc.list(
        lapply(fittedModel$runjags$mcmc, corSamples,
               covData=cdat,
               thetaUnique=fittedModel$mptInfo$thetaUnique,
               rho=ifelse(class(fittedModel) == "betaMPT", TRUE, FALSE),
               corProbit = fittedModel$mptInfo$corProbit))
    }
  }

  fittedModel$mcmc.summ <- summarizeMCMC(tmp)
  fittedModel$summary <- summarizeMPT(mcmc = tmp,
                                      summ = fittedModel$mcmc.summ,
                                      mptInfo = fittedModel$mptInfo)
  fittedModel$call <- c(fittedModel$call, match.call())
  fittedModel
}

