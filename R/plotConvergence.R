
#### wrappers for convenient convergence plots

#' @export
#' @describeIn plot Plot convergence for beta MPT
plot.betaMPT <- function(x, parameter="mean", type="default", ...){
  plot.traitMPT(x,parameter=parameter, type=type,...)
}


#' Plot Convergence for Hierarchical MPT Models
#'
#' @param x fitted hierarchical MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param parameter which parameter to plot (e.g., \code{"theta"}, \code{"rho"}, \code{"slope"}). Parameters are matched partially, in order to plot all entries of vector valued parameters (see \code{\link{getParam}} to get a list of parameters)
#' @param type type of convergence plot. Can be one of \code{"default"} (trace+density), \code{"trace"}, \code{"autocorr"}, \code{"crosscorr"},\code{"density"},   \code{"gelman"}. See, e.g., \code{\link[coda]{plot.mcmc.list}}
#' @param ... further arguments passed to the plotting functions in coda
#' @export
#' @describeIn plot Plot convergence for latent-trait MPT
#' @importFrom coda traceplot acfplot gelman.plot as.mcmc.list varnames crosscorr.plot autocorr.plot densplot
plot.traitMPT <- function(x, parameter="mean", type="default", ...){

  #   if(type %in% c("ecdf", "histogram", "autocorr", "key", "crosscorr","all")){
  #     plot(x$fit, plot.type = type, vars = parameter,...)
  #   }else{

  MPT.mcmc <- as.mcmc.list(x$mcmc) # as.mcmc.list(x$mcmc$BUGSoutput)
  allnam <- varnames(MPT.mcmc)
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
  coda::varnames(MPT.mcmc) <- allnam
  # `varnames()<-`(MPT.mcmc, allnam)

  switch(type,
         "trace" = traceplot(MPT.mcmc[,idx],...),
         "acf" = acfplot(MPT.mcmc[,idx],...),
         "gelman" = gelman.plot(MPT.mcmc[,idx],...),
         "crosscorr" = crosscorr.plot(MPT.mcmc[,idx],...),
         "autocorr" = autocorr.plot(MPT.mcmc[,idx],...),
         "density" = densplot(MPT.mcmc[,idx],...),
         "default" = plot(MPT.mcmc[,idx],...),
         stop("Check 'type' for possible plots." )
  )

}
