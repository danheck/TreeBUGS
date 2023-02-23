
#' Plot Prior Distributions
#'
#' Plots prior distributions for group means, standard deviation, and correlations of MPT parameters across participants.
#' @param M number of random samples to approximate priors of group-level parameters
#' @param probitInverse which latent-probit parameters (for \code{\link{traitMPT}} model) to transform to probability scale. Either \code{"none"}, \code{"mean"} (simple transformation \eqn{\Phi(\mu)}), or \code{"mean_sd"} (see \code{\link{probitInverse}})
#' @param ... further arguments passed to \code{plot}
#' @inheritParams priorPredictive
#' @details This function samples from a set of hyperpriors (either for hierarchical traitMPT or betaMPT structure) to approximate the implied prior distributions on the parameters of interest (group-level mean, SD, and correlations of MPT parameters). Note that the normal distribution \code{"dnorm(mu,prec)"} is parameterized as in JAGS by the mean and precision (= 1/variance).
#' @export
#' @seealso \code{\link{priorPredictive}}
#' @examples
#' \dontrun{
#' # default priors for traitMPT:
#' plotPrior(list(
#'   mu = "dnorm(0,1)", xi = "dunif(0,10)",
#'   V = diag(2), df = 2 + 1
#' ), M = 4000)
#'
#' # default priors for betaMPT:
#' plotPrior(list(
#'   alpha = "dgamma(1,.1)",
#'   beta = "dgamma(1,.1)"
#' ), M = 4000)
#' }
plotPrior <- function(prior, probitInverse = "mean", M = 5000, nCPU = 3, ...) {
  ############### prior samples
  samples <- sampleHyperprior(prior, M, # S=1,
    probitInverse = probitInverse, truncSig = .995, nCPU = nCPU
  )
  model <- attr(samples, "model")
  S <- ncol(samples$mean)
  if (model == "traitMPT") {
    S.plot <- ifelse(S > 1 && (max(length(prior$xi), length(prior$mu)) > 1 |
      any(prior$V != diag(S))), S, 1)
  } else {
    S.plot <- ifelse(S > 1 && (max(length(prior$alpha), length(prior$beta)) > 1), S, 1)
  }

  ################# plotting

  mfrow <- par()$mfrow
  qq <- seq(0, 1, .05)
  bins <- min(60, round(M / 40))
  histcol <- adjustcolor("gray", alpha.f = .7)
  par(mfrow = c(2, ifelse(model == "traitMPT", 2, 1)))

  ######################## MEAN
  for (s in 1:S.plot) {
    hist(samples$mean[, s], bins,
      freq = FALSE, col = histcol,
      main = paste0(
        "Prior on group mean: ",
        ifelse(model == "traitMPT",
          paste0("mu=", ifelse(length(prior$mu) == 1, prior$mu[1], prior$mu[s])),
          paste0("alpha=", ifelse(length(prior$alpha) == 1, prior$alpha[1], prior$alpha[s]))
        )
      ),
      xlab = paste0(
        "Group mean",
        ifelse(model == "traitMPT" && probitInverse == "none",
          " (probit scale)",
          " (probability scale)"
        )
      ),
      border = histcol, las = 1, ...
    )
  }

  ######################## SD
  for (s in 1:S.plot) {
    hist(samples$sd[, s], bins,
      freq = FALSE, col = histcol,
      xlim = c(0, max(max(samples$sd, na.rm = TRUE), .5)),
      main = paste0(
        "Prior on group SD: ",
        ifelse(model == "traitMPT",
          paste0("xi=", ifelse(length(prior$xi) == 1, prior$xi[1], prior$xi[s])),
          paste0("beta=", ifelse(length(prior$beta) == 1, prior$beta[1], prior$beta[s]))
        )
      ),
      xlab = paste0(
        "Group SD",
        ifelse(model == "betaMPT" || probitInverse == "mean_sd",
          " (probability scale)",
          " (probit scale)"
        )
      ),
      border = histcol, las = 1, ...
    )
  }

  ######################## CORRELATION
  if (model == "traitMPT" && S > 1) {
    for (s1 in 1:(S - 1)) {
      for (s2 in (s1 + 1):S) {
        hist(samples$rho[s1, s2, ], bins,
          freq = FALSE, col = histcol,
          main = paste0("Correlation (", s1, " and ", s2, "): df=", prior$df),
          xlab = "Correlation (probit scale)", border = histcol, las = 1, ...
        )
      }
    }
  }

  par(mfrow = mfrow)
}
