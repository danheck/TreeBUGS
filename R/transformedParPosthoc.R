
#' Get Transformed Parameters
#'
#' Computes transformations of MPT parameters based on the MCMC posterior samples (e.g., differences of parameters)
#' @inheritParams posteriorPredictive
#' @param transformedParameters list with parameter transformations that should be computed based on the posterior samples (e.g., for testing parameter differences: \code{list("diffD=Do-Dn")})
#' @param level whether to compute transformations of \code{"group"} or \code{"individual"} estimates
#' @param nCPU number of CPU cores across which the MCMC chains are distributed
#' @return an \link[coda]{mcmc.list} of posterior samples for the transformed parameters
#' @examples
#' \dontrun{
#' tt <- transformedParameters(fittedModel,
#'                             list("diff=a-b","p=a>b"),
#'                             level="individual")
#' summary(tt)
#' }
#' @export
transformedParameters <- function (fittedModel,
                                   transformedParameters,
                                   level = "group",
                                   nCPU = 4){

  parsed <- getTransformed(fittedModel$mptInfo$thetaNames,
                           transformedParameters=transformedParameters,
                           mergeString=FALSE)
  S <- length(fittedModel$mptInfo$thetaUnique)
  P <- length(parsed$transformedParameters)
  N <- nrow(fittedModel$mptInfo$data)

  samp <- list()
  #### for a single MCMC chain:
  getTrans <- function(mm){
    if(level == "group"){
      samp <- mm[,rep(1,P)]
      colnames(samp) <- parsed$transformedParameters
      for(i in 1:P){
        samp[,i] <- apply(mm[,paste0("mean[",1:S,"]")], 1,
                                 function(xx) eval(parse(text=parsed$modelstring[1]),
                                                   list(mean=unlist(xx))))
      }
    }else{
      samp <- mm[,rep(1,P*N)]
      for(i in 1:P){
        idx <- (i-1)*N+1:N
        colnames(samp)[idx] <- paste0(parsed$transformedParameters[i],"[", 1:N,"]")
        for(n in 1:N){
          samp[,idx[n]] <- apply(mm[,paste0("theta[",1:S,",",n,"]")], 1,
                                        function(xx) eval(parse(text=parsed$modelstring[1]),
                                                          list(mean=unlist(xx))))
        }
      }
    }
    samp
  }

  cl <- makeCluster(nCPU)
  clusterExport(cl, c("S","N","P","parsed",
                      "transformedParameters","level"), envir=environment())
  # loop across MCMC chains:
  samp <- parLapply(cl, fittedModel$runjags$mcmc, getTrans)
  stopCluster(cl)

  as.mcmc.list(samp)
}

