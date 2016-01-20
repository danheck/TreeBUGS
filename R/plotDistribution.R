
#' Plot Distribution of Individual Estimates
#'
#' Plots histograms of the mean posterior estimates for participants including the estimated shape on the group level (e.g., the beta distribution in case of the beta-MPT)
#'
#' @inheritParams plotFit
#' @param scale only for latent-trait MPT: should estimates be plotted on the \code{"latent"} or the \code{"probability"} scale (i.e., as MPT parameters). Can be abbreviated by \code{"l"}  and \code{"p"}.
#' @param ... further arguments passed to \code{\link{hist}} (e.g., \code{breaks=50} to get a more fine-grained histogram or \code{xlim=0:1} to use the same x-scale for all histograms)
#' @seealso \code{\link{plot.traitMPT}}
#' @export
plotDistribution <- function(fittedModel, scale="probability", ...){
  mfrow <- par()$mfrow
  mar <- par()$mar
  scale <- match.arg(scale, c("probability", "latent"))

  means <- fittedModel$summary$groupParameters$mean[,1]
  parnames <- names(fittedModel$summary$individParameters[,1,1])

  S <- length(means)
  nrow <- floor(sqrt(S))
  ncol <- ceiling(sqrt(S))
  par(mfrow=c(min(4, nrow), min(6, ncol)), mar=c(2, 2, 3, .3))

  for(idx in 1:S ){
    indEsts <- fittedModel$summary$individParameters[idx,,1]
    if(class(fittedModel) == "traitMPT"){
      sigma <- fittedModel$mcmc$BUGSoutput$mean$sigma

      if(scale == "latent"){
        hist(qnorm(indEsts), freq=F, main=parnames[idx], xlab="Latent scale", ...)
        xx <- seq(0, 1, .005)
        lines(xx, dnorm(xx, qnorm(means[idx]), sigma[idx]), col=2)
      }else{
        hist(indEsts, freq=F, main=parnames[idx], xlab="Probability scale", ...)
        # values on latent scale:
        xx <- seq(-5, 5, .005)
        # values on probability scale:
        xx.p <- pnorm(xx)
        # discrete approximation to density on latent scale:
        p.diff <- diff(c(pnorm(xx, qnorm(means[idx]), sigma[idx]), 1))
        lines(xx.p, p.diff/diff(c(xx.p,1)), col=2)
      }

    }else if(class(fittedModel) == "betaMPT"){
      alpha <- fittedModel$summary$groupParameters$alpha[,1]
      beta <- fittedModel$summary$groupParameters$beta[,1]

      hist(indEsts, freq=F, main=parnames[idx], xlab="Probability scale", ...)
      xx <- seq(0, 1, .005)
      lines(xx, dbeta(xx, alpha[idx], beta[idx]), col=2)

    }
  }

  par(mfrow=mfrow, mar=mar)
}
