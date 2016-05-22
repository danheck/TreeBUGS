

#' Compute Posterior Predictive P-Values
#'
#' Computes posterior predictive p-values to test model fit.
#' @inheritParams posteriorPredictive
#' @author Daniel Heck
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#' @export
PPP <- function(fittedModel, M=1000, nCPU=4){

  cats <- fittedModel$mptInfo$MPT$Category
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  numItems <-   t(apply(fittedModel$mptInfo$data, 1,
                        function(x)  tapply(x, fittedModel$mptInfo$MPT$Tree, sum)))
  # selection list for mapping of categories to trees:
  sel <- lapply(TreeNames, function(tt) tree %in% tt)

  # expected frequencies:
  freq.exp <- posteriorPredictive(fittedModel, M, expected = TRUE, nCPU=nCPU)
  M <- length(freq.exp)

  cl <- makeCluster(nCPU)
  clusterExport(cl, c("TreeNames"), envir=environment())
  # sample conditional on expected probabilities:
  freq.pred <- lapply(freq.exp, function(fe){
    for(k in 1:length(TreeNames)){
      fe[,sel[[k]]] <- t(apply(fe[,sel[[k]]], 1,
                               function(x) rmultinom(1, size=sum(x), prob=x/sum(x))))
    }
    fe
  })
  stopCluster(cl)

  freq.obs <- fittedModel$mptInfo$data[,colnames(freq.pred[[1]])]

  # mean frequencies:
  mean.obs <- colMeans(freq.obs)
  mean.pred <- t(sapply(freq.pred, colMeans))
  mean.exp <- t(sapply(freq.exp, colMeans))

  # statistics:
  T1.obs <- T1.pred <-T2.obs <- T2.pred <- NULL
  try({
    T1.obs <- apply(mean.exp, 1, T1stat, n=mean.obs)
    T2.obs <- sapply(freq.exp, T2stat, n.ind=freq.obs, tree=tree)
    for(m in 1:M){
      T1.pred[m] <- T1stat(mean.exp[m,], mean.pred[m,])
      T2.pred[m] <- T2stat(n.ind.exp=freq.exp[[m]],
                           n.ind=freq.pred[[m]],
                           tree=tree)
    }
  })

  # PPP-value:
  T1.p <- mean(T1.obs < T1.pred)
  T2.p <- mean(T2.obs < T2.pred)
  res <- list(T1.obs=T1.obs, T1.pred=T1.pred,
              T2.obs=T2.obs, T2.pred=T2.pred,
              T1.p=T1.p, T2.p=T2.p)
  class(res) <- "postPredP"
  res
}

#' @export
print.postPredP <- function(x, ...){
  cat(" ## Mean structure (T1):\n",
      "Observed = ", mean(x$T1.obs), "; Predicted = ", mean(x$T1.pred), "; p-value = ", mean(x$T1.p),"\n",
      "## Covariance structure (T2):\n",
      "Observed = ", mean(x$T2.obs), "; Predicted = ", mean(x$T2.pred), "; p-value = ", mean(x$T2.p),"\n")
}

T1stat <- function(n.exp, n){
  n.exp[n.exp==0] <- 1e-10
  dev <- (n-n.exp)^2/n.exp
  sum(dev)
}



# check:
# n.ind.exp <- freq.exp[[1]]
# n.ind <- freq.obs

#' @importFrom  stats cov
#' @importFrom  utils combn
T2stat <- function(n.ind.exp, n.ind, tree){

  N <- nrow(n.ind)
  ncat <- ncol(n.ind)

  # expected (sigma)
  sigma <- cov(n.ind.exp)
  for(k in 1:length(unique(tree))){
    sel <- tree %in% unique(tree)[k]
    tmp <- n.ind.exp[,sel]
    n <- rowSums(tmp)  #*(N-1)/N

    # variance for multinomial: n*p1*(1-p1)
    mean.ind.cov <-  diag(colMeans(tmp*(1-tmp/n)  ))

    # covariance for multinomial: -n*p1*p2
    combs <- combn(1:sum(sel),2)
    for(c in 1:ncol(combs)){
      mean.ind.cov[combs[1,c],combs[2,c]] <- mean(-tmp[,combs[1,c]]*tmp[,combs[2,c]]/n)
    }
    mean.ind.cov[lower.tri(mean.ind.cov)] <- mean.ind.cov[upper.tri(mean.ind.cov)]

    # add to covariance of expected frequencies:
    sigma[sel,sel] <-  sigma[sel,sel] + (N-1)/N*mean.ind.cov
  }


  # observed
  sd <- cov(n.ind)

  # discrepancy:
  dev <- (sd - sigma)^2
  denom <-  1/sqrt(diag(sigma))
  normalize <- matrix(denom %x% denom, ncat)
  sum(dev*normalize)
}
