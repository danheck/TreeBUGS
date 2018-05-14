
############## PLOTTING FUNCTIONS FOR betaMPT

#' Plot Parameter Estimates
#'
#' Plot parameter estimates for hierarchical MPT models.
#'
#' @param x a fitted Beta or latent-trait MPT model
#' @param includeIndividual whether to plot individual estimates
#' @param addLines whether to connect individual parameter estimates by lines
#' @param estimate type of point estimates for group-level and individual parameters
#'     (either \code{"mean"} or \code{"median"})
#' @param select character vector of parameters to be plotted (e.g., \code{select = c("d", "g")}. Can be used to plot subsets of parameters and change the order of parameters.
#' @param ... further arguments passed to the standard \code{\link{plot}} function
#'
#' @author Daniel Heck
#' @seealso \code{\link{betaMPT}}, \code{\link{traitMPT}}, \code{\link{plotDistribution}}
#' @examples
#' \dontrun{
#' plotParam(fit, addLines = TRUE,
#'           estimate = "median",
#'           select = c("d1", "d2"))
#' }
#' @export
plotParam <- function(x, includeIndividual = TRUE, addLines = FALSE,
                      estimate = "mean", select = "all", ...){

  stat <- ifelse(estimate == "median", "50%", "Mean")
  par.group <- x$summary$groupParameters$mean
  par.ind <- x$summary$individParameters
  parnames <- substr(rownames(par.group), 6, 100)
  if (select[1] == "all"){
    select <- parnames
  } else {
    if (!all(select %in% parnames))
      stop("Check arguments: Not all parameters in 'select' are included in the MPT model!\n",
           "Parameters are: ", paste(parnames, collapse=", "))
    par.group <- par.group[paste0("mean_", select),, drop = FALSE]
  }
  dims <- dim(par.ind)
  S <- nrow(par.group)  # parameters
  N <- dims[2]          # persons
  means <- par.group[,stat]

  par(mfrow=c(1,1))
  plot(1:S, means, ylim=0:1, xlim=c(.5, S+.5),pch=19, xaxt="n", #size=3,
       xlab = "MPT Parameters",
       ylab=paste0("Estimate (",estimate ,"s)"), col=2,
       main=paste0("Group-level ",estimate,"s + 95% CI (red)",
                   ifelse(includeIndividual,
                          paste0(" and individual ",estimate,"s (gray)"),
                          "")), ...)
  axis(side = 1, at = 1:S, labels=select)
  if(includeIndividual){
    for(i in 1:N){
      if (addLines){
        lines(1:S + .05, par.ind[select,i,stat],
              col = adjustcolor(col = "black", alpha.f = .5))
        points(1:S + .05, par.ind[select,i,stat], cex=.9, pch=16,
              col = adjustcolor(col = "black", alpha.f = .5))
      } else {
        points(1:S + seq(-.2,.2, length.out = N)[i],
               col = adjustcolor(col = "black", alpha.f = .5), #col=rainbow(N, alpha=.4)[i],
               pch = 16,
               par.ind[select,i,stat], cex=.9)
      }
    }
    points(1:S, means, cex=1.3,  col=2,pch=19)
  }
  segments(x0=1:S, y0=par.group[,3],
           y1=par.group[,5], lwd=2, col=2)

}
