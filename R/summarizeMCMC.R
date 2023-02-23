#' MCMC Summary
#'
#' TreeBUGS-specific MCMC summary for \code{mcmc.list}-objects.
#'
#' @param mcmc a \code{\link[coda]{mcmc.list}} object
#' @param batchSize size of batches of parameters used to reduce memory load
#'   when computing posterior summary statistics (including Rhat and effective
#'   sample size).
#' @param probs quantile probabilities used to compute credibility intervals
#'
#' @importFrom coda varnames
#' @export
summarizeMCMC <- function(
    mcmc,
    batchSize = 50,
    probs = c(.025, .50, .975)
) {
  if (inherits(mcmc, c("traitMPT", "betaMPT", "simpleMPT"))) {
    mcmc <- mcmc$runjags$mcmc
  }
  if (inherits(mcmc, "runjags")) {
    mcmc <- mcmc$mcmc
  }

  # initialize matrix with summary statistics
  vnames <- varnames(mcmc)
  npar <- length(vnames)
  snames <- c(
    "Mean", "SD", names(quantile(1, probs)),
    "Time-series SE", "n.eff", "Rhat", "R_95%"
  )
  summTab <- matrix(NA, npar, length(snames),
    dimnames = list(vnames, snames)
  )

  # summarize in batches (to avoid RAM issues)
  n_batches <- (npar %/% batchSize) + 1
  for (ii in seq(n_batches)) {
    if (n_batches == 1) {
      idx <- seq(npar)
    } else if (ii < n_batches) {
      idx <- (ii - 1) * batchSize + 1:batchSize # complete batches
    } else if ((ii - 1) * batchSize + 1 <= npar) {
      idx <- ((ii - 1) * batchSize + 1):npar
    } else {
      break()
    }

    try({
      mcmc.mat <- do.call("rbind", mcmc[, idx, drop = FALSE])
      summTab[idx, "Mean"] <- apply(mcmc.mat, 2, mean, na.rm = TRUE)
      summTab[idx, "SD"] <- apply(mcmc.mat, 2, sd, na.rm = TRUE)
      summTab[idx, 2 + seq(length(probs))] <-
        t(apply(mcmc.mat, 2, quantile, probs, na.rm = TRUE))
      rm(mcmc.mat)
      gc(verbose = FALSE)
    })

    try(
      {
        summTab[idx, "n.eff"] <- round(effectiveSize(mcmc[, idx]))
        summTab[idx, "Time-series SE"] <- summTab[idx, "SD"] / sqrt(summTab[idx, "n.eff"])
      },
      silent = TRUE
    )

    try(summTab[idx, c("Rhat", "R_95%")] <-
      gelman.diag(mcmc[, idx], multivariate = FALSE)[[1]])
  }

  if (all(is.na(summTab))) {
    cat("summarizeMCMC: posterior summary in baches failed. trying coda::summary instead.\n")
    try({
      summ <- summary(mcmc)
      summTab <- cbind(summ[[1]][, c("Mean", "SD")],
        summ[[2]],
        "Time-series SE" = summ[[1]][, "Time-series SE"],
        "n.eff" = (summ[[1]][, "SD"] / summ[[1]][, "Time-series SE"])^2,
        "Rhat" = NA, "R_95%" = NA
      )
      ###### MCMC effective N:
      # "Time-series SE" = "SD" / sqrt("Effective N")
      # "Effective N"    = ("SD" / "Time-series SE")^2
    })
  }

  summTab
}
