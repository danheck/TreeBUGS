
############## PLOTTING FUNCTIONS FOR betaMPT

#' Plot Parameter Estimates
#'
#' Plot parameter estimates for hierarchical MPT models
#'
#' @param x a fitted Beta or latent-trait MPT model
#' @param includeIndividual whether to plot individual estimates
#' @param estimate type of point estimates for group-level and individual parameters (either \code{"mean"} or \code{"median"})
#' @param ... further arguments passed to the standard \code{\link{plot}} function
#' @author Daniel Heck
#' @seealso \code{\link{betaMPT}}, \code{\link{traitMPT}}, \code{\link{plotDistribution}}
#' @export
plotParam <- function(x, includeIndividual=TRUE, estimate = "mean", ...){

  stat <- ifelse(estimate=="median","50%","Mean")
  dims <- dim(x$summary$individParameters)
  N <- dims[2]
  S <- dims[1]
  means <- x$summary$groupParameters$mean[,stat]
  means <- x$summary$groupParameters$mean[,"50%"]

  par(mfrow=c(1,1))
  plot(1:S, means, ylim=0:1, xlim=c(.5, S+.5),pch=19, xaxt="n", #size=3,
       xlab = "MPT Parameters", ylab="Estimate", col=2,
       main=paste0("Group-level ",estimate,"s (+ 95% CI)",
                   ifelse(includeIndividual,paste0(" and individual ",estimate,"s"),"")), ...)
  axis(side = 1, at = 1:S, labels=substr(names(means), 6, 100))
  segments(x0=1:S, y0=x$summary$groupParameters$mean[,3],
           y1=x$summary$groupParameters$mean[,5], lwd=2, col=2)
  if(includeIndividual){
    for(i in 1:N){
      points(1:S+seq(-.2,.2, length.out = N)[i],
             col=adjustcolor(col = "black", alpha.f=.5), #col=rainbow(N, alpha=.4)[i],
             pch=16,
             x$summary$individParameters[,i,stat], cex=.9)
    }
    points(1:S, means, cex=1.3,  col=2,pch=19)
  }

}
