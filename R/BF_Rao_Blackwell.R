#' Bayes Factors for Simple (Nonhierarchical) MPT Models
#'
#' Computes Bayes factors for simple (fixed-effects, nonhierarchical) MPT models
#' with beta distributions as priors on the parameters.
#'
#' @param models list of models fitted with \code{\link{simpleMPT}}, e.g.,
#'   \code{list(mod1, mod2)}
#' @param resample how many of the posterior samples of the MPT parameters
#'   should be resampled per model
#' @param store whether to save parameter samples
#' @inheritParams marginalMPT
#'
#' @details Currently, this is only implemented for a single data set!
#'
#' Uses a Rao-Blackwellized version of the product-space method (Carlin & Chib,
#' 1995) as proposed by Barker and Link (2013). First, posterior distributions
#' of the MPT parameters are approximated by independent beta distributions.
#' Second, for one a selected model, parameters are sampled from these proposal
#' distributions. Third, the conditional probabilities to switch to a different
#' model are computed and stored. Finally, the eigenvector with eigenvalue one
#' of the matrix of switching probabilities provides an estimate of the
#' posterior model probabilities.
#'
#' @seealso \code{\link{marginalMPT}}
#' @references Barker, R. J., & Link, W. A. (2013). Bayesian multimodel
#'   inference by RJMCMC: A Gibbs sampling approach. The American Statistician,
#'   67(3), 150-156.
#'
#'   Carlin, B. P., & Chib, S. (1995). Bayesian model choice via Markov chain
#'   Monte Carlo methods. Journal of the Royal Statistical Society. Series B
#'   (Methodological), 57(3), 473-484.
#' @export
BayesFactorMPT <- function(
    models,
    dataset = 1,
    resample,
    batches = 5,
    scale = 1,
    store = FALSE,
    cores = 1
) {
  if (!is.list(models) || any(sapply(models, class) != "simpleMPT")) {
    stop("'models' must be a list models with fitted simpleMPT!")
  }
  datas <- lapply(models, function(m) m$mptInfo$data)
  M <- length(models)
  if (M < 2) {
    stop("At least two models must be provided.")
  }
  for (i in 2:M) {
    if (any(datas[[1]] != datas[[i]])) {
      stop("each model must have one vector of frequencies that must be identical for all ")
    }
  }
  if (missing(resample)) {
    resample <- min(sapply(models, function(x) {
      length(x$runjags$mcmc) * nrow(x$runjags$mcmc[[1]])
    }))
  }

  # 2. Approximate posteriors by beta densities
  betapars <- shape.prior <- list()
  for (m in 1:M) {
    ab <- approximatePosterior(models[[m]], dataset = dataset, sample = 500)
    betapars[[m]] <- pmax(scale * ab, .1)

    prior_pars <- do.call("cbind", models[[m]]$mpt$hyperprior)
    shape.prior[[m]] <- prior_pars[rownames(betapars[[m]]), , drop = FALSE]
  }

  # 3. Loop 1 (rows of P): Model k
  #     => sample palette vector: phi = c(t1, ..., tk, ..., tm)
  #     => tk ~ resample from MCMC posteriors
  #     => ti ~ sample from approximation

  P <- array(NA, c(M, M, resample))

  # m = row of P (FROM which model to start jumping)
  row.P <- function(m) {
    # theta <- vector("list", M)
    # theta[[m]] <- resampling(models[[m]], dataset = dataset, resample=resample)
    # theta[-m] <- lapply(betapars[-m],
    #                     function(bp) apply(bp, 1,
    #                                        function(s) rbeta(resample, s[1], s[2])))
    # loglik <- mapply(llMPT, pars = theta, mod = models, MoreArgs = list(dataset = dataset))
    # prior.pseudo <- mapply(dProductBeta, x = theta, shapes = betapars)
    # prior.current <- mapply(dProductBeta, x = theta, shapes = shape.prior)
    loglik <- prior.pseudo <- prior.current <- matrix(NA, resample, M)
    for (i in 1:M) {
      if (i == m) {
        theta <- resampling(models[[m]], dataset = dataset, resample = resample)
      } else {
        theta <- rProductBeta(resample, betapars[[i]])
      }
      loglik[, i] <- llMPT(
        pars = theta, mod = models[[i]],
        dataset = dataset
      )
      prior.pseudo[, i] <- dProductBeta(x = theta, shapes = betapars[[i]])
      prior.current[, i] <- dProductBeta(x = theta, shapes = shape.prior[[i]])
    }

    ### 4. Loop 2 (entries in row k): Compute transition probabilities
    ###     => P(i|k) = P(y|tk,Mk) * P(tk|Mk) * prod(P(ti|Mk)) * P(Mk)
    # prior for full palette vector:
    posterior <-
      exp(loglik + # y | theta_k, M_k
        prior.current + # theta_k | M_k
        rowSums(prior.pseudo) - prior.pseudo) # theta_i | M_k
    posterior / rowSums(posterior)
  }
  if (cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, c("resample", "M", "models", "betapars", "shape.prior", "dataset"),
      envir = environment()
    )
    P.tmp <- parSapply(cl, 1:M, row.P, simplify = FALSE)
    stopCluster(cl)
  } else {
    P.tmp <- sapply(1:M, row.P, simplify = FALSE)
  }
  for (m in 1:M) {
    P[m, , ] <- t(P.tmp[[m]])
  }
  rm(P.tmp) # ; gc()

  #### 5. Rao-Blackwell Estimate: Average probabilities
  ####    Left eigenvector with eigenvalue = 1
  P.mean <- apply(P, 1:2, mean)
  ev2 <- rb.estimate(P.mean)
  # Re(eigen(t(P.mean))$vec[,1])
  p.est <- rb.estimate(P.mean) # ev2/sum(ev2)

  # batch estimate + SE
  idx <- rep(1:batches, each = round(resample / batches))
  if (length(idx) > resample) idx <- idx[1:resample]
  if (length(idx) < resample) idx <- c(idx, rep(batches, resample - length(idx)))
  tmp <- matrix(NA, batches, M)
  for (b in 1:batches) {
    P.mean <- apply(P[, , idx == b], 1:2, mean)
    # ev2 <- Re(eigen(t(P.mean))$vec[,1])
    tmp[b, ] <- rb.estimate(P.mean) # ev2/sum(ev2)
  }
  batchSE <- apply(tmp, 2, sd) / sqrt(batches)
  p.batch <- matrix(c(colMeans(tmp), batchSE),
    nrow = 2, ncol = M, byrow = TRUE,
    dimnames = list(c("Mean", "SE"), Model = names(models))
  )

  dimnames(P.mean) <- list("from" = names(models), "to" = names(models))
  names(p.est) <- names(models)
  res <- list(
    "posterior" = p.est,
    "p.batch" = p.batch,
    "P.mean" = P.mean
  )
  if (store) {
    res$samples <- P
  }
  res
}

rb.estimate <- function(P) {
  ev <- Re(eigen(t(P))$vec[, 1])
  ev / sum(ev)
}
