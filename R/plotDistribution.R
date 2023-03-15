#' Plot Distribution of Individual Estimates
#'
#' Plots histograms of the posterior-means of individual MPT parameters against
#' the group-level distribution given by the posterior-mean of the hierarchical
#' parameters (e.g., the beta distribution in case of the beta-MPT)
#'
#' @param fittedModel fitted latent-trait or beta MPT model
#'   (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param scale only for latent-trait MPT: should estimates be plotted on the
#'   \code{"latent"} or the \code{"probability"} scale (i.e., as MPT
#'   parameters). Can be abbreviated by \code{"l"}  and \code{"p"}.
#' @param ... further arguments passed to \code{\link{hist}} (e.g.,
#'   \code{breaks=50} to get a more fine-grained histogram)
#'
#' @details For the latent-trait MPT, differences due to continuous predictors
#'   or discrete factors are currently not considered in the group-level
#'   predictions (red density). Under such a model, individual estimates are not
#'   predicted to be normally distributed on the latent scale as shown in the
#'   plot.
#'
#' @seealso \code{\link{plot.traitMPT}}
#' @export
plotDistribution <- function(
    fittedModel,
    scale = "probability",
    ...
) {
  mfrow <- par()$mfrow
  mar <- par()$mar
  scale <- match.arg(scale, c("probability", "latent"))

  means <- fittedModel$summary$groupParameters$mean[, 1]
  parnames <- names(fittedModel$summary$individParameters[, 1, 1])

  S <- length(means)
  nrow <- floor(sqrt(S))
  ncol <- ceiling(sqrt(S))
  par(mfrow = c(min(4, nrow), min(6, ncol)), mar = c(2, 2, 3, .3))

  for (idx in 1:S) {
    indEsts <- fittedModel$summary$individParameters[idx, , 1]
    if (all(is.na(indEsts)))
      stop("No MCMC samples for the individual-level MPT parameters (theta) were stored. \n",
           "Please re-fit model with the argument:  monitorIndividual = TRUE")

    if (inherits(fittedModel, "traitMPT")) {
      sigma <- fittedModel$summary$groupParameters$sigma[, "Mean"]
      # sigma <- fittedModel$mcmc$BUGSoutput$mean$sigma

      # values on latent scale:
      xx <- seq(-10, 10, length.out = 3000)
      if (scale == "latent") {
        hist(qnorm(indEsts),
          freq = F, main = paste0("Parameter ", parnames[idx]),
          col = "gray", xlab = "Latent scale", las = 1, ...
        )
        lines(xx, dnorm(xx, qnorm(means[idx]), sigma[idx]), col = 2)
      } else {
        hist(indEsts,
          freq = F, main = paste0("Parameter ", parnames[idx]), xlim = 0:1,
          col = "gray", xlab = "Probability scale", las = 1, ...
        )
        # values on probability scale:
        xx.p <- pnorm(xx)
        # discrete approximation to density on latent scale:
        p.diff <- diff(c(pnorm(xx, qnorm(means[idx]), sigma[idx]), 1))
        lines(xx.p, p.diff / diff(c(xx.p, 1)), col = 2)
      }
    } else if (inherits(fittedModel, "betaMPT")) {
      alpha <- fittedModel$summary$groupParameters$alpha[, 1]
      beta <- fittedModel$summary$groupParameters$beta[, 1]

      hist(indEsts,
        freq = F, main = paste0("Parameter ", parnames[idx]), xlim = 0:1,
        col = "gray", xlab = "Probability scale", las = 1, ...
      )
      xx <- seq(0, 1, length.out = 1000)
      lines(xx, dbeta(xx, alpha[idx], beta[idx]), col = 2)
    }
  }

  par(mfrow = mfrow, mar = mar)
}
