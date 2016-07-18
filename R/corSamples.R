
# rho: whether to compute theta inter-correlations (for beta MPT)
# corProbit: whether to use theta parameters on probit scale
corSamples <- function(mcmc, covData, thetaUnique,
                       rho=FALSE, corProbit = FALSE){

  # extract theta samples:
  S <- length(thetaUnique)
  sel.theta <- setdiff(grep("theta", varnames(mcmc), fixed=TRUE),
                       grep("thetaFE", varnames(mcmc), fixed=TRUE))
  theta <- mcmc[,sel.theta, drop = FALSE]
  if(corProbit){
    theta <- qnorm(theta)
    tmin <- min(theta[theta != -Inf])
    tmax <- max(theta[theta != Inf])
    theta[theta == -Inf] <- min(-5, tmin)
    theta[theta == Inf] <- max(5, tmax)
  }

  ########## correlations: theta <-> covData
  if(!is.null(covData) && !ncol(covData) == 0){
    cor.samp <- matrix(t(apply(theta, 1, function(tt){
      theta.tmp <- matrix(tt, nrow=nrow(covData), byrow=TRUE)
      cor(theta.tmp, covData)
    })), nrow = nrow(theta))
    cor.names <- outer(thetaUnique, colnames(covData), FUN=paste, sep="_")
    colnames(cor.samp) <- paste0("cor_", cor.names)
  }else{
    cor.samp <- NULL
  }

  ########## intercorrelations for theta (=rho)
  if(rho & ncol(theta)>1){
    rho.samp <- matrix(t(apply(theta, 1, function(tt){
      theta.tmp <- matrix(tt, nrow=length(tt)/S, byrow=TRUE)
      cc <- cor(theta.tmp)
      # cc[lower.tri(cc)]
      cc
    })), nrow = nrow(theta))
    rho.names <- outer(1:length(thetaUnique),
                       1:length(thetaUnique), FUN=paste, sep=",")
    # colnames(rho.samp) <- paste0("rho[",rho.names[lower.tri(rho.names)],"]")
    colnames(rho.samp) <- paste0("rho[",rho.names,"]")
  }else{
    rho.samp <- NULL
  }

  mcmc.attr <- attr(mcmc, "mcpar")
  mcmc(cbind(mcmc, cor.samp, rho.samp),
       start=mcmc.attr[1], end = mcmc.attr[2], thin=mcmc.attr[3])
}
