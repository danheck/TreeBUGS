

# T1 statistic:


getPPP <- function(fittedModel, M=1000, nCPU=4){

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
  T1.obs <- T1.pred <- apply(mean.exp, 1, T1stat, n=mean.obs)
  for(m in 1:M){
    T1.pred[m] <- T1stat(mean.exp[m,], mean.pred[m,])
  }

  # PPP-value:
  T1.p <- mean(T1.obs < T1.pred)
  list(T1.obs=T1.obs, T1.pred=T1.pred, T1.p=T1.p)
}


T1stat <- function(n.exp, n){
  dev <- (n-n.exp)^2/n.exp
  sum(dev)
}

