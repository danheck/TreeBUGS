
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
#' @param parameter which parameter to plot (e.g., \code{"theta"}, \code{"mean"}, \code{"rho"}, \code{"slope"}). Parameters are matched partially, in order to plot all entries of vector valued parameters (see \code{\link{getParam}} to get a list of parameters). Moreover, parameter labels can be used, e.g., \code{"theta[D]"} or \code{"rho[D,g]"}
#' @param type type of convergence plot. Can be one of \code{"default"} (trace+density), \code{"acf"} (auto-correlation function), \code{"trace"}, \code{"autocorr"}, \code{"crosscorr"},\code{"density"},   \code{"gelman"}. See, e.g., \code{\link[coda]{mcmc.list}}
#' @param ... further arguments passed to the plotting functions in coda
#' @export
#' @describeIn plot Plot convergence for latent-trait MPT
#' @importFrom coda traceplot acfplot gelman.plot as.mcmc.list varnames crosscorr.plot autocorr.plot densplot
plot.traitMPT <- function(x, parameter="mean", type="default", ...){
  mcmc <- x$runjags$mcmc
  allnam <- varnames(mcmc)
  thetaUnique <- x$mptInfo$thetaUnique
  S <- length(thetaUnique)
  parameter <- gsub(" ", "", parameter, fixed = TRUE)

  # unnecessary rho-parameters
  if (parameter == "rho"){
    rho.idx <- outer(1:S, 1:S, paste, sep = ",")
    rho.double <- paste0("rho[", rho.idx[lower.tri(rho.idx, diag = TRUE)], "]")
    allnam <- setdiff(allnam, rho.double)
    mcmc <- mcmc[,allnam]
  }
  parameter <- name2idx(parameter, thetaUnique)

  idx <- setdiff(grep(parameter, allnam, fixed = TRUE),
                 grep(".pred", allnam, fixed = TRUE))
  if(length(idx) <=0){
    stop("Parameter not found in MCMC object.")
  }
  if(parameter == "theta"){
    allnam[idx] <-  paste0(allnam[idx], rep(thetaUnique,"_", length(idx)/2))
  } else {
    allnam <- idx2name(allnam, thetaUnique)
  }
  coda::varnames(mcmc) <- allnam

  switch(type,
         "trace" = traceplot(mcmc[,idx],...),
         "acf" = acfplot(mcmc[,idx],...),
         "gelman" = gelman.plot(mcmc[,idx],...),
         "crosscorr" = crosscorr.plot(mcmc[,idx],...),
         "autocorr" = autocorr.plot(mcmc[,idx],...),
         "density" = densplot(mcmc[,idx],...),
         "default" = plot(mcmc[,idx],...),
         stop("Check 'type' for possible plots." )
  )

}

idx2name <- function(parnames, thetaUnique){
  for (i in seq_along(thetaUnique)){
    parnames <- gsub(paste0("[", i),
                     paste0("[", thetaUnique[i]),
                     parnames, fixed = TRUE)
    parnames <- gsub(paste0(i, "]"),
                     paste0(thetaUnique[i], "]"),
                     parnames, fixed = TRUE)
  }
  parnames
}

name2idx <- function(parnames, thetaUnique){
  for (i in seq_along(thetaUnique)){
    parnames <- gsub(paste0("[", thetaUnique[i]),
                     paste0("[", i),
                     parnames, fixed = TRUE)
    parnames <- gsub(paste0(thetaUnique[i], "]"),
                     paste0(i, "]"),
                     parnames, fixed = TRUE)
  }
  parnames
}
