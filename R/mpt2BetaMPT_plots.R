
############## PLOTTING FUNCTIONS FOR betaMPT
#' Plot Beta-MPT Parameter Estimates
#'
#' @param x a fitted Beta-MPT model (see \link{mpt2BetaMPT})
#' @param includeIndividual whether to plot individual estimates
#' @param ... further arguments passed to the standard \link{plot} function
#' @export
plot.betaMPT <- function(x, includeIndividual=TRUE, ...){

  if(! x$sampler %in% c("JAGS", "jags")){
    warning("Adjusted plotting functions for hierarchical MPT models only available when using JAGS.")
  }else{
    dims <- dim(x$summary$individParameters)
    N <- dims[2]
    S <- dims[1]
    means <- x$summary$meanParameters$mean[,1]

    plot(1:S, means, ylim=0:1, pch=19, xaxt="n", size=1.5,
         xlab = "MPT Parameters", ylab="Estimate",
         main="Mean Estimates (Including 95% Credibility Interval of Estimate)", ...)
    axis(side = 1, at = 1:S, labels=substr(names(means), 6, 100))
    segments(x0=1:S, y0=x$summary$meanParameters$mean[,3],
             y1=x$summary$meanParameters$mean[,5], lwd=1.5)
    if(includeIndividual){
      for(i in 1:N){
        points(1:S+runif(1,-.1,.1), col=rainbow(N, alpha=.4)[i], pch=16,
               x$summary$individParameters[,i,1], size=.9)
      }
    }
  }


}
