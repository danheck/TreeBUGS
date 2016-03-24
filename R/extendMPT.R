
#' Extend MCMC Sampling for MPT Model
#'
#' @param fittedModel a fitted \code{\link{traitMPT}} or \code{\link{betaMPT}}
#' @inheritParams betaMPT
#' @param ... further arguments passed to \link[runjags]{extend.jags}
#' @export
extendMPT <- function(fittedModel, n.iter = 10000, n.adapt = 1000, n.burnin = 0, n.thin = 1, ...){

  tmp <- extend.jags(fittedModel$runjags,
                     burnin = n.burnin,
                     sample = ceiling((n.iter-n.burnin)/n.thin),
                     adapt = n.adapt,
                     thin=n.thin,summarise = FALSE, ...)

  fittedModel$mcmc.summ <- summarizeMCMC(tmp$mcmc)
  fittedModel$summary <- summarizeMPT(mcmc = tmp$mcmc,
                                      summ = fittedModel$mcmc.summ,
                                      mptInfo = fittedModel$mptInfo)
  fittedModel$call <- c(fittedModel$call, match.call())
  fittedModel
}
