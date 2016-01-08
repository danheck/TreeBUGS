
############## PLOTTING FUNCTIONS FOR betaMPT
#' Plot Beta-MPT Parameter Estimates
#'
#' @param x a fitted Beta-MPT model (see \code{\link{betaMPT}})
#' @param includeIndividual whether to plot individual estimates
#' @param ... further arguments passed to the standard \code{\link{plot}} function
#' @author Daniel Heck
#' @export
plot.betaMPT <- function(x, includeIndividual=TRUE, ...){

  if(! x$sampler %in% c("JAGS", "jags")){
    warning("Adjusted plotting functions for hierarchical MPT models only available when using JAGS.")
  }else{
    dims <- dim(x$summary$individParameters)
    N <- dims[2]
    S <- dims[1]
    means <- x$summary$groupParameters$mean[,1]

    plot(1:S, means, ylim=0:1, xlim=c(.5, S+.5),pch=19, xaxt="n", #size=3,
         xlab = "MPT Parameters", ylab="Estimate",
         main="Mean estimates (including 95% credibility interval for group mean)", ...)
    axis(side = 1, at = 1:S, labels=substr(names(means), 6, 100))
    segments(x0=1:S, y0=x$summary$groupParameters$mean[,3],
             y1=x$summary$groupParameters$mean[,5], lwd=2)
    if(includeIndividual){
      for(i in 1:N){
        points(1:S+seq(-.2,.2, length.out = N)[i], col=rainbow(N, alpha=.4)[i], pch=16,
               x$summary$individParameters[,i,1], cex=.9)
      }
      points(1:S, means, cex=1.3,pch=19)
    }
  }


}

#' Plots for latent-trait MPT models (see \code{\link{traitMPT}}).
#'
#' @inheritParams  plot.betaMPT
#' @author Daniel Heck
#' @export
plot.traitMPT <- function(x, includeIndividual=TRUE, ...){
  plot.betaMPT(x,...)
}


# convergence diagnostic: pass to R2JAGS?

# #' @export
# traceplot.betaMPT <- function(x,...){traceplot(x$mcmc,...)}
#
# #' @export
# autocorr.plot.betaMPT <- function(x,...){autocorr.plot(x$mcmc,...)}
