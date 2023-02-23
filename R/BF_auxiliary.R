
# Approximate posterior distribution of MPT parameters by a well-known, simple density
# Curently, only distribution="beta" is implemented
approximatePosterior <- function(mod, dataset = 1,
                                 sample = 500, distribution = "beta",
                                 lower = .1, upper = 1e4) {
  # estimate alpha/beta parameters of beta approximation
  S <- length(mod$mptInfo$thetaUnique)
  betapar <- matrix(1, S, 2, dimnames = list(
    mod$mptInfo$thetaUnique,
    c("alpha", "beta")
  ))
  for (i in 1:S) {
    sel <- paste0("theta[", i, ",", dataset, "]")
    ss <- unlist(mod$runjags$mcmc[, sel])
    ss <- sample(ss, min(sample, length(ss)))
    m <- mean(ss)
    v <- var(ss)
    betapar[i, 1] <- m * (m * (1 - m) / v - 1)
    betapar[i, 2] <- (1 - m) * (m * (1 - m) / v - 1)
    try(
      betapar[i, ] <- fitdistr(ss, "beta",
        list(
          shape1 = betapar[i, 1],
          shape2 = betapar[i, 2]
        ),
        lower = lower,
        upper = upper
      )$estimate,
      silent = TRUE
    )
  }
  betapar
}


# resample MCMC iterations from simpleMPT object
resampling <- function(mod, dataset = 1, resample = 1000) {
  S <- length(mod$mptInfo$thetaUnique)

  sel <- paste0("theta[", 1:S, ",", dataset, "]")
  C <- length(mod$runjags$mcmc)
  R <- nrow(mod$runjags$mcmc[[1]])
  r <- ceiling(resample / C)
  if (resample > C * R) {
    warning(
      "Fitted models have less samples than required for resampling.",
      "Posterior samples will be reused!"
    )
    rr <- lapply(
      mod$runjags$mcmc[, sel, drop = FALSE],
      function(mm) mm[sample(1:R, r, replace = TRUE), , drop = FALSE]
    )
  } else {
    rr <- lapply(
      mod$runjags$mcmc[, sel, drop = FALSE],
      function(mm) mm[sample(1:R, r), , drop = FALSE]
    )
  }
  do.call("rbind", rr)[1:resample, , drop = FALSE]
}
