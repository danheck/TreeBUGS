
#' C++ Sampler for Hierarchical Beta-MPT Model
#'
#' Fast Gibbs sampler in C++ that is tailored to the beta-MPT model.
#'
#' @inheritParams betaMPT
#' @param shape shape parameter of Gamma-hyperdistribution for \eqn{\alpha_s} and \eqn{\beta_s}
#' @param rate rate parameter of Gamma-hyperdistribution
#' @author Daniel Heck
#' @importFrom parallel parLapply stopCluster detectCores
#' @export
betaMPTcpp <- function(eqnfile, data, restrictions,
                       covData, corProbit=FALSE,
                       n.iter=20000, n.burnin = 2000,
                       n.thin = 5,  n.chains=3, ppp = 0,
                       shape = 1, rate = .1, parEstFile){

  if(missing(restrictions)) restrictions <- NULL
  if(missing(covData)) covData <- NULL

  tab <- readEQN(eqnfile)
  mpt <- parseEQN(tab)

  mpt.res <- parseRestrictions(mpt, restrictions)
  thetaUnique <- colnames(mpt.res$a)

  data <- readData(data)
  data <- data[,mpt.res$cat.names]
  N <- nrow(data)
  covData <- covDataRead(covData, N)

  time0 <- Sys.time()
  cat("MCMC sampling started at ", format(time0), "\n")

  mcmc.list <- list()

  sim <- function(idx){
    sim <- betampt(M = n.iter, H = as.matrix(data),
                   a = mpt.res$a, b = mpt.res$b, c = mpt.res$c, map = mpt.res$map,
                   shape = shape, rate=rate )

    S <- ncol(mpt.res$a)
    if(S>1){
      colnames(sim$mean) <- paste0("mean[",1:S,"]")
      colnames(sim$sd) <- paste0("sd[",1:S,"]")
      colnames(sim$alph) <- paste0("alph[",1:S,"]")
      colnames(sim$bet) <- paste0("bet[",1:S,"]")
    }
    tnames <- outer(1:S,paste0(",",1:N), paste0)
    tt <- matrix(sim$theta, nrow = n.iter, ncol=S*N,
                 dimnames=list(NULL, paste0("theta[",c(t(tnames)),"]")))
    tmp <- with(sim, cbind(mean, sd, alph, bet, tt ))
    if(S == 1)
      colnames(tmp)[1:4] <- c("mean","sd","alph","bet")
    mcmc(tmp, start = n.burnin+1, end=n.iter, thin = n.thin)
  }

  ncpu <- min(detectCores(), n.chains)
  cl <- makeCluster(ncpu)
  # tttt <- clusterEvalQ(cl, library(TreeBUGS))
  clusterExport(cl, c("sim","n.iter","data","mpt.res",
                      "shape","rate","N","n.burnin","n.thin"), envir = environment())
  mcmc.list <- parLapply(cl, 1:n.chains, sim)
  stopCluster(cl)
  # for(cc in 1:n.chains){
  #   sim <- betampt(M = n.iter, H = as.matrix(data),
  #                  a = mpt.res$a, b = mpt.res$b, c = mpt.res$c, map = mpt.res$map,
  #                  shape = shape, rate=rate )
  #
  #   S <- ncol(mpt.res$a)
  #   if(S>1){
  #     colnames(sim$mean) <- paste0("mean[",1:S,"]")
  #     colnames(sim$sd) <- paste0("sd[",1:S,"]")
  #     colnames(sim$alph) <- paste0("alph[",1:S,"]")
  #     colnames(sim$bet) <- paste0("bet[",1:S,"]")
  #   }
  #   tnames <- outer(1:S,paste0(",",1:N), paste0)
  #   tt <- matrix(sim$theta, nrow = n.iter, ncol=S*N,
  #                dimnames=list(NULL, paste0("theta[",c(t(tnames)),"]")))
  #   tmp <- with(sim, cbind(mean, sd, alph, bet, tt ))
  #   if(S == 1)
  #     colnames(tmp)[1:4] <- c("mean","sd","alph","bet")
  #   mcmc.list[[cc]] <- mcmc(tmp, start = n.burnin+1, end=n.iter, thin = n.thin)
  # }
  mcmc.list <- as.mcmc.list(mcmc.list)
  time1 <- Sys.time()
  cat("MCMC sampling finished at", format(time1), "\n  ")
  print(time1-time0)

  # store details about model:
  hyperprior <- list(alpha = check.hyperprior(paste0("dgamma(",shape,",",rate,")"), thetaUnique, label="alpha"),
                     beta = check.hyperprior(paste0("dgamma(",shape,",",rate,")"), thetaUnique, label="beta"))
  mptInfo <- list(model="betaMPT",
                  thetaNames = data.frame(Parameter = colnames(mpt.res$a),
                                          theta=1:ncol(mpt.res$a)),
                  thetaUnique = thetaUnique,
                  thetaFixed = NULL,
                  MPT=mpt.res,
                  eqnfile=eqnfile,
                  data=data,
                  restrictions=restrictions,
                  covData=covData,
                  corProbit=corProbit,
                  predTable=NULL,
                  predFactorLevels=NULL,
                  predType=NULL,
                  transformedParameters=NULL, #transformedPar$transformedParameters,
                  hyperprior=hyperprior)

  # correlation of posterior samples:
  if(!is.null(covData)){
    sel <-  is.numeric( data.frame(covData)[1,])
    cdat <- as.matrix(covData[,sel,drop = FALSE])
  }else{
    cdat  <- NULL
  }
  mcmc.list <- as.mcmc.list(
    lapply(mcmc.list, corSamples,
           covData=cdat, thetaUnique=thetaUnique,
           rho=TRUE, corProbit = corProbit))



  # own summary (more stable than runjags)
  mcmc.summ <- summarizeMCMC(mcmc.list)
  summary <- summarizeMPT(mcmc = mcmc.list, mptInfo = mptInfo,
                          summ=mcmc.summ)
  summary$dic <- NULL

  # class structure for TreeBUGS
  fittedModel <- list(summary=summary,
                      mptInfo=mptInfo,
                      mcmc.summ = mcmc.summ,
                      runjags = list(mcmc=mcmc.list),
                      postpred=NULL,
                      call=match.call(),
                      time=time1-time0)
  class(fittedModel) <- "betaMPT"

  if(ppp>0){
    cat("\nComputing posterior-predictive p-values....\n")
    postPred <- PPP(fittedModel, M=ppp,
                    nCPU=length(mcmc.list))
    fittedModel$postpred <- postPred[c("freq.exp", "freq.pred", "freq.obs")]
    try(
      fittedModel$summary$fitStatistics <- list(
        "overall"=c(
          "T1.observed"=mean(postPred$T1.obs),
          "T1.predicted"=mean(postPred$T1.pred),
          "p.T1"=postPred$T1.p,
          "T2.observed"=mean(postPred$T2.obs),
          "T2.predicted"=mean(postPred$T2.pred),
          "p.T2"=postPred$T2.p
        ))
    )
  }

  # write results to file
  writeSummary(fittedModel, parEstFile)

  # varnames(m2$runjags$mcmc)
  ## TODO: Rcpp parallel

  fittedModel
}
