
#### wrappers for convenient convergence plots

#' @export
#' @describeIn plot Plot convergence for beta MPT
plot.betaMPT <- function(x, parameter="mean", type="default", ...){
  plot.traitMPT(x,parameter=parameter, type=type,...)
}

#' @export
#' @describeIn plot Plot convergence for nonhierarchical MPT model
plot.simpleMPT <- function(x, type="default", ...){
  plot.traitMPT(x,parameter="theta", type=type,...)
}


#' Plot Convergence for Hierarchical MPT Models
#'
#' @param x fitted hierarchical MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param parameter which parameter to plot (e.g., \code{"theta"}, \code{"rho"}, \code{"slope"}). Parameters are matched partially, in order to plot all entries of vector valued parameters (see \code{\link{getParam}} to get a list of parameters)
#' @param type type of convergence plot. Can be one of \code{"default"} (trace+density), \code{"acf"} (auto-correlation function), \code{"trace"}, \code{"autocorr"}, \code{"crosscorr"},\code{"density"},   \code{"gelman"}. See, e.g., \code{\link[coda]{plot.mcmc.list}}
#' @param ... further arguments passed to the plotting functions in coda
#' @export
#' @describeIn plot Plot convergence for latent-trait MPT
#' @importFrom coda traceplot acfplot gelman.plot as.mcmc.list varnames crosscorr.plot autocorr.plot densplot
plot.traitMPT <- function(x, parameter="mean", type="default", ...){

  #   if(type %in% c("ecdf", "histogram", "autocorr", "key", "crosscorr","all")){
  #     plot(x$fit, plot.type = type, vars = parameter,...)
  #   }else{

  # MPT.mcmc <- x$runjags$mcmc
  allnam <- varnames(x$runjags$mcmc)
  thetaUnique <- x$mptInfo$thetaUnique
  idx <- setdiff(grep(parameter,allnam) , grep(".pred",allnam))
  if(length(idx) <=0){
    stop("Parameter not found in MCMC object.")
  }
  # parLabel  <- names(x$mcmc$BUGSoutput$mean[[parameter]]) # names(x$mcmc$BUGSoutput$mean[[parameter]])
  if(length(idx) == length(thetaUnique)){
    allnam[idx] <- paste0(parameter, "_", thetaUnique)
  }else if(parameter == "theta"){
    allnam[idx] <-  paste0(allnam[idx], rep(thetaUnique,"_", length(idx)/2))
  }
  coda::varnames(x$runjags$mcmc) <- allnam
  # `varnames()<-`(MPT.mcmc, allnam)

  switch(type,
         "trace" = traceplot(x$runjags$mcmc[,idx],...),
         "acf" = acfplot(x$runjags$mcmc[,idx],...),
         "gelman" = gelman.plot(x$runjags$mcmc[,idx],...),
         "crosscorr" = crosscorr.plot(x$runjags$mcmc[,idx],...),
         "autocorr" = autocorr.plot(x$runjags$mcmc[,idx],...),
         "density" = densplot(x$runjags$mcmc[,idx],...),
         "default" = plot(x$runjags$mcmc[,idx],...),
         stop("Check 'type' for possible plots." )
  )

}
