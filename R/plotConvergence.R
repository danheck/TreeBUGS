
#### wrappers for convenient convergence plots

#' @export
#' @describeIn plot Plot convergence for beta MPT
plot.betaMPT <- function(x, parameter="mean", type="trace", ...){
  plot.traitMPT(x,...)
}


#' Plot Convergence for Hierarchical MPT Models
#'
#' @param x fitted hierarchical MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param parameter which parameter to plot (e.g., \code{"theta"}, \code{"rho"}, \code{"slope"}). Parameters are matched partially, in order to plot all entries of vector valued parameters (see \code{\link{getParam}} to get a list of parameters)
#' @param type what type of convergence plot (\code{"trace"}, \code{"acf"}, \code{"gelman"}). See, e.g., \code{\link[coda]{traceplot}}
#' @param ... further arguments passed to the \code{coda}-plotting function
#' @export
#' @describeIn plot Plot convergence for latent-trait MPT
#' @importFrom coda traceplot acfplot gelman.plot as.mcmc.list varnames
plot.traitMPT <- function(x, parameter="mean", type="trace", ...){

  MPT.mcmc <- as.mcmc.list(x$mcmc$BUGSoutput)
  allnam <- varnames(MPT.mcmc)
  idx <- setdiff(grep(parameter,allnam) , grep(c("response."),allnam))
  if(length(idx) <=0){
    stop("Parameter not found in MCMC object.")
  }
  parLabel  <- names(x$mcmc$BUGSoutput$mean[[parameter]])
  allnam[idx] <- paste0(parameter, "_", parLabel)
  coda::varnames(MPT.mcmc) <- allnam
  # `varnames()<-`(MPT.mcmc, allnam)

  switch(type,
         "trace" = traceplot(MPT.mcmc[,idx],...),
         "acf" = acfplot(MPT.mcmc[,idx],...),
         "gelman" = gelman.plot(MPT.mcmc[,idx],...))
}
