#' Compute Posterior Predictive P-Values
#'
#' Computes posterior predictive p-values to test model fit.
#'
#' @inheritParams posteriorPredictive
#' @param T2 whether to compute T2 statistic to check coveriance structure (can
#'   take a lot of time). If some participants do not have responses for some
#'   trees, (co)variances are computed by pairwise deletion of the corresponding
#'   persons.
#' @param type whether the T1 statistic of expected means is computed using
#'   Person's \code{"X2"} or the likelihood-ratio statistic \code{"G2"}
#'
#' @author Daniel Heck
#' @references Klauer, K. C. (2010). Hierarchical multinomial processing tree
#' models: A latent-trait approach. Psychometrika, 75, 70-98.
#'
#' @importFrom  stats cov
#' @importFrom  utils combn
#' @importFrom parallel clusterMap
#' @export
PPP <- function(
    fittedModel,
    M = 1000,
    nCPU = 4,
    T2 = TRUE,
    type = "X2"
){

  cats <- fittedModel$mptInfo$MPT$Category
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  numItems <- t(apply(
    fittedModel$mptInfo$data, 1,
    function(x) tapply(x, tree, sum)
  ))
  # selection list for mapping of categories to trees:
  sel <- lapply(TreeNames, function(tt) tree %in% tt)

  # expected frequencies:
  freq.exp <- posteriorPredictive(fittedModel, M, expected = TRUE, nCPU = nCPU)
  M <- length(freq.exp)




  # sample conditional on expected probabilities: (TODO: Cpp sampling; cf. mpt2irt)
  freq.pred <- freq.exp
  for (r in 1:length(freq.pred)) {
    for (k in 1:length(TreeNames)) {
      for (i in 1:nrow(freq.pred[[r]])) {
        nn <- freq.exp[[r]][i, sel[[k]]]
        if (sum(nn) > 0) {
          freq.pred[[r]][i, sel[[k]]] <- rmultinom(1, round(sum(nn)), nn)
        } else {
          freq.pred[[r]][i, sel[[k]]] <- 0
        }
      }
    }
  }

  freq.obs <- fittedModel$mptInfo$data[, colnames(freq.pred[[1]])]

  # mean frequencies:
  mean.obs <- colMeans(freq.obs)
  mean.pred <- t(sapply(freq.pred, colMeans))
  mean.exp <- t(sapply(freq.exp, colMeans))

  cl <- makeCluster(nCPU)
  clusterExport(cl, c(
    "TreeNames", "T1stat", "mean.obs", "mean.pred", "mean.exp",
    "freq.exp", "freq.obs", "freq.pred", "type"
  ),
  envir = environment()
  )

  # statistics:
  T1.obs <- T1.pred <- T2.obs <- T2.pred <- NULL
  ind.T1.obs <- ind.T1.pred <- matrix(NA, M, nrow(freq.obs))
  try(
    {
      T1.obs <- parApply(cl, mean.exp, 1, T1stat, n = mean.obs, type = type)
      # T1stat(mean.exp[10,], mean.obs) ; T1.obs[10]
      T1.pred <- clusterMap(cl, T1stat,
        n.exp = matAsList(mean.exp),
        n = matAsList(mean.pred),
        MoreArgs = list(type = type),
        SIMPLIFY = TRUE
      )
      # T1stat(mean.exp[1,], mean.pred[1,]) ; T1.pred[1]

      # individual T1:
      ind.T1.pred <- do.call("rbind", mapply(T1stat.mat,
        n.exp = freq.exp, n = freq.pred,
        MoreArgs = list(type = type), SIMPLIFY = FALSE
      ))
      ind.T1.obs <- do.call("rbind", mapply(T1stat.mat,
        n.exp = freq.exp, n = list(freq.obs)[rep(1, M)],
        MoreArgs = list(type = type), SIMPLIFY = FALSE
      ))
      # ind.T1.obs <- t(sapply(freq.exp, T1stat, n = list(freq.obs)[rep(1,M)] ))
      # mm <- 15; i <- 16
      # T1stat(n.exp = freq.exp[[mm]][i,], n = freq.pred[[mm]][i,]) ; ind.T1.pred[mm,i]
      # T1stat(n.exp = freq.exp[[mm]][i,], n = freq.obs[i,]) ; ind.T1.obs[mm,i]

      if (T2) {
        T2.obs <- parSapply(cl, freq.exp, T2stat, n.ind = freq.obs, tree = tree)
        # T2stat(freq.exp[[1]], n.ind=freq.obs, tree) ; T2.obs[1]
        T2.pred <- clusterMap(cl, T2stat,
          n.ind.exp = freq.exp,
          n.ind = freq.pred,
          MoreArgs = list(tree = tree), SIMPLIFY = TRUE
        )
        # T2stat(freq.exp[[1]], n.ind=freq.pred[[1]], tree); T2.pred[1]
      }
    },
    silent = TRUE
  )

  stopCluster(cl)


  # PPP-value:
  T1.p <- mean(T1.obs < T1.pred)
  T2.p <- mean(T2.obs < T2.pred)
  ind.T1.p <- colMeans(ind.T1.obs < ind.T1.pred)
  names(ind.T1.p) <- 1:nrow(freq.obs)
  res <- list(
    T1.obs = T1.obs, T1.pred = T1.pred, T1.p = T1.p,
    T2.obs = T2.obs, T2.pred = T2.pred, T2.p = T2.p,
    ind.T1.obs = ind.T1.obs,
    ind.T1.pred = ind.T1.pred,
    ind.T1.p = ind.T1.p,
    freq.exp = freq.exp, freq.pred = freq.pred, freq.obs = freq.obs
  )
  class(res) <- "postPredP"
  res
}

addPPP <- function(fittedModel, M = 1000) {
  if (M > 0) {
    cat("\nComputing posterior-predictive p-values....\n")
    postPred <- PPP(fittedModel,
      M = M,
      nCPU = length(fittedModel$runjags$mcmc)
    )
    fittedModel$postpred <- postPred[c("freq.exp", "freq.pred", "freq.obs")]
    try(
      fittedModel$summary$fitStatistics <- list(
        "overall" = c(
          "T1.observed" = mean(postPred$T1.obs),
          "T1.predicted" = mean(postPred$T1.pred),
          "p.T1" = postPred$T1.p,
          "T2.observed" = mean(postPred$T2.obs),
          "T2.predicted" = mean(postPred$T2.pred),
          "p.T2" = postPred$T2.p
        ),
        "individual" = list(
          "T1.obs" = postPred$ind.T1.obs,
          "T1.pred" = postPred$ind.T1.pred,
          "T1.p" = postPred$ind.T1.p
        )
      )
    )
  }
  fittedModel
}

#' @export
print.postPredP <- function(x, ...) {
  cat(
    " ## Mean structure (T1):\n",
    "Observed = ", mean(x$T1.obs),
    "; Predicted = ", mean(x$T1.pred),
    "; p-value = ", mean(x$T1.p), "\n"
  )
  if (!is.null(x$T2.obs)) {
    cat(
      "\n## Covariance structure (T2):\n",
      "Observed = ", mean(x$T2.obs),
      "; Predicted = ", mean(x$T2.pred),
      "; p-value = ", mean(x$T2.p), "\n"
    )
  }

  cat("\n## Individual fit (T1):\n")
  print(round(x$ind.T1.p, 3))
}

matAsList <- function(matrix) {
  lapply(apply(matrix, 1, list), function(x) x[[1]])
}


T1stat.mat <- function(n.exp, n, type = "X2", const = 1e-8) {
  n.exp[n.exp == 0] <- const
  if (type == "X2") {
    dev <- (t(n) - t(n.exp))^2 / t(n.exp)
  } else if (type == "G2") {
    dev <- 2 * t(n) * log(t(n) / t(n.exp))
    dev[t(n) == 0] <- 0
  }
  colSums(dev)
}

T1stat <- function(n.exp, n, type = "X2", const = 1e-8) {
  n.exp[n.exp == 0] <- const
  if (type == "X2") {
    dev <- (n - n.exp)^2 / n.exp
  } else if (type == "G2") {
    dev <- 2 * n * log(n / n.exp)
    dev[n == 0] <- 0
  }
  sum(dev)
}

T2stat <- function(n.ind.exp, n.ind, tree) {
  N <- nrow(n.ind)
  ncat <- ncol(n.ind)

  # expected (sigma)
  sigma <- cov(n.ind.exp)
  for (k in 1:length(unique(tree))) {
    sel <- tree %in% unique(tree)[k]
    tmp <- n.ind.exp[, sel]
    n <- rowSums(tmp) #* (N-1)/N

    # variance for multinomial: n*p1*(1-p1)
    varmat <- tmp * (1 - tmp / n)
    # varmat[tmp == 0] <- 0
    mean.ind.cov <- diag(colMeans(varmat, na.rm = TRUE))

    # covariance for multinomial: -n*p1*p2
    combs <- combn(1:sum(sel), 2)
    for (c in 1:ncol(combs)) {
      mean.ind.cov[combs[1, c], combs[2, c]] <- mean(-tmp[, combs[1, c]] * tmp[, combs[2, c]] / n, na.rm = TRUE)
    }
    mean.ind.cov[lower.tri(mean.ind.cov)] <- t(mean.ind.cov)[lower.tri(mean.ind.cov)]

    # add to covariance of expected frequencies:
    sigma[sel, sel] <- sigma[sel, sel] + (N - 1) / N * mean.ind.cov
  }


  # observed
  sd <- cov(n.ind)

  # discrepancy:
  dev <- (sd - sigma)^2
  denom <- 1 / sqrt(diag(sigma))
  normalize <- matrix(denom %x% denom, ncat)
  sum(dev * normalize)
}
