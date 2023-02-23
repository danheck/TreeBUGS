
#' Marginal Likelihood for Simple MPT
#'
#' Computes the marginal likelihood for simple (fixed-effects, nonhierarchical) MPT models.
#'
#' @inheritParams simpleMPT
#' @param method either \code{"importance"} (importance sampling using a mixture of uniform and beta-aproximation of the posterior) or \code{"prior"} (brute force Monte Carlo sampling from prior)
#' @param posterior number of posterior samples used to approximate importance-sampling densities (i.e., beta distributions)
#' @param dataset for which data set should Bayes factors be computed?
#' @param mix mixture proportion of the uniform distribution for the importance-sampling density
#' @param scale how much should posterior-beta approximations be downscaled to get fatter importance-sampling density
#' @param samples total number of samples from parameter space
#' @param batches number of batches. Used to compute a standard error of the estimate.
#' @param show whether to show progress
#' @param cores number of CPUs used
#'
#' @details
#' Currently, this is only implemented for a single data set!
#'
#' If \code{method = "prior"}, a brute-force Monte Carlo method is used and parameters are directly sampled from the prior.Then, the likelihood is evaluated for these samples and averaged (fast, but inefficient).
#'
#' Alternatively, an importance sampler is used if \code{method = "importance"}, and the posterior distributions of the MPT parameters are approximated by independent beta distributions. Then each parameter \eqn{s} is sampled from the importance density:
#'
#' \eqn{mix*U(0,1) + (1-mix)*Beta(scale*a_s, scale*b_s)}
#'
#' @examples
#' # 2-High-Threshold Model
#' eqn <- "## 2HTM ##
#'    Target  Hit  d
#'    Target  Hit  (1-d)*g
#'    Target  Miss (1-d)*(1-g)
#'    Lure    FA   (1-d)*g
#'    Lure    CR   (1-d)*(1-g)
#'    Lure    CR   d"
#' data <- c(
#'   Hit = 46, Miss = 14,
#'   FA = 14, CR = 46
#' )
#'
#' # weakly informative prior for guessing
#' aa <- c(d = 1, g = 2)
#' bb <- c(d = 1, g = 2)
#' curve(dbeta(x, aa["g"], bb["g"]))
#'
#' # compute marginal likelihood
#' htm <- marginalMPT(eqn, data,
#'   alpha = aa, beta = bb,
#'   posterior = 200, samples = 1000
#' )
#' # second model: g=.50
#' htm.g50 <- marginalMPT(eqn, data, list("g=.5"),
#'   alpha = aa, beta = bb,
#'   posterior = 200, samples = 1000
#' )
#'
#' # Bayes factor
#' # (per batch to get estimation error)
#' bf <- htm.g50$p.per.batch / htm$p.per.batch
#' mean(bf) # BF
#' sd(bf) / sqrt(length(bf)) # standard error of BF estimate
#'
#' @seealso \code{\link{BayesFactorMPT}}
#' @references
#' Vandekerckhove, J. S., Matzke, D., & Wagenmakers, E. (2015). Model comparison and the principle of parsimony. In Oxford Handbook of Computational and Mathematical Psychology (pp. 300-319). New York, NY: Oxford University Press.
#' @export
marginalMPT <- function(eqnfile, data, restrictions, alpha = 1, beta = 1,
                        dataset = 1, method = "importance", posterior = 500,
                        mix = .05, scale = .9, samples = 10000, batches = 10,
                        show = TRUE, cores = 1) {
  if (mix < 0 | mix > 1 | scale < 0 | scale > 1) {
    stop("The tuning parameters 'mix' and 'scale' must be in the interval [0,1].")
  }
  if (is.character(data)) {
    data <- readData(data)
  }
  if (!is.vector(data) && nrow(data) > 1) {
    data <- data[dataset, ]
  }

  ############### 1. fit MPT / get MPT structure
  tmp <- capture.output(
    mod <- simpleMPT(
      eqnfile = eqnfile, data = data, restrictions = restrictions,
      n.iter = ifelse(method == "importance", posterior + 200, 5),
      n.burnin = ifelse(method == "importance", 200, 2),
      n.thin = 1, n.chains = 1, alpha = alpha, beta = beta
    )
  )

  if (method == "importance") {
    betapar <- approximatePosterior(mod) * scale
    sampling <- sample.importance
  } else if (method == "prior") {
    sampling <- sample.prior
    betapar <- NULL
  }
  t0 <- Sys.time()
  if (show) cat("Sampling parameters to estimate marginal likelihood...\n")

  n.per.batch <- ceiling(samples / batches)
  if (cores > 1) {
    cl <- makeCluster(cores)
    samp <- parSapply(cl, 1:batches, sampling,
      mod = mod, betapar = betapar,
      mix = mix, samples = n.per.batch, show = show
    )
    stopCluster(cl)
  } else {
    samp <- sapply(1:batches, sampling,
      mod = mod, betapar = betapar,
      mix = mix, samples = n.per.batch, show = show
    )
  }
  t1 <- Sys.time()

  ############### 4. Aproximate integral + batch SE
  p <- c(mean = mean(c(samp)), SE = sd(c(samp)) / sqrt(length(samp)))
  p.batch <- colMeans(samp)
  p.batch <- c(
    "p" = mean(p.batch),
    "SE_p" = sd(p.batch) / sqrt(batches),
    "log_p" = mean(log(p.batch)),
    "SE_log_p" = sd(log(p.batch)) / sqrt(batches)
  )

  res <- list(
    "p" = p,
    "p.batch" = p.batch,
    "p.per.batch" = colMeans(samp),
    "prior" = do.call("cbind", mod$mptInfo$hyperprior),
    "data" = mod$mptInfo$data,
    "sampler" = list(
      "method" = method,
      "samples" = samples,
      "time" = t1 - t0
    )
  )
  if (method == "importance") {
    res$sampler$mix <- mix
    res$sampler$scale <- scale
  }
  class(res) <- "marginalMPT"
  if (show) {
    cat("\nFinished in ", format(t1 - t0))
  }
  return(res)
}


# mixture: mix*U(0,1) + (1-mix)*Beta(a,b)
# sampling:
rmix <- function(i, samples, betapar, mix) {
  n.unif <- rbinom(1, samples, mix)
  c(
    runif(n.unif),
    rbeta(samples - n.unif, betapar[i, 1], betapar[i, 2])
  )
}
# density:
dmix <- function(x, betapar, mix) {
  sum(log(mix + (1 - mix) * dbeta(x, betapar[, 1], betapar[, 2])))
}


sample.importance <- function(b, samples, mod,
                              betapar, mix, show = FALSE) {
  S <- length(mod$mptInfo$thetaUnique)
  const <- logMultinomCoefficient(mod)
  if (show) cat(b, "")
  xx <- sapply(1:nrow(betapar), rmix,
    samples = samples, betapar = betapar, mix = mix
  )
  gx <- apply(xx, 1, dmix, betapar = betapar, mix = mix)
  # prior:
  px <- rep(0, samples)
  for (s in 1:S) {
    px <- px + dbeta(xx[, s],
      log = TRUE,
      shape1 = mod$mptInfo$hyperprior$alpha[s],
      shape2 = mod$mptInfo$hyperprior$beta[s]
    )
  }
  # likelihood:
  fx <- c(loglikMPT(xx,
    h = unlist(mod$mptInfo$data),
    a = mod$mptInfo$MPT$a,
    b = mod$mptInfo$MPT$b,
    c = mod$mptInfo$MPT$c,
    map = mod$mptInfo$MPT$map
  ))
  # importance weights wx <- px-gx
  exp(fx + px - gx + const)
}

sample.prior <- function(b, samples, mod,
                         betapar = NULL, mix = NULL, show = FALSE) {
  S <- length(mod$mptInfo$thetaUnique)
  const <- logMultinomCoefficient(mod)
  if (show) cat(b, "")
  xx <- matrix(NA, samples, S)
  for (s in 1:S) {
    xx[, s] <- rbeta(samples,
      shape1 = mod$mptInfo$hyperprior$alpha[s],
      shape2 = mod$mptInfo$hyperprior$beta[s]
    )
  }
  # likelihood:
  fx <- c(loglikMPT(xx,
    h = unlist(mod$mptInfo$data),
    a = mod$mptInfo$MPT$a,
    b = mod$mptInfo$MPT$b,
    c = mod$mptInfo$MPT$c,
    map = mod$mptInfo$MPT$map
  ))
  exp(fx + const)
}
