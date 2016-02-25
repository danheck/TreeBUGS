
#' Get Mean Group Estimates
#'
#' For hierarchical latent-trait MPT models with discrete predictor variables.
#' @param traitMPT a fitted \code{\link{traitMPT}} model
#' @param factor whether to get group estimates for all combinations of factor levels (default) or only for specific factors (requires the names of the covariates in covData)
# @param computeT1 whether to compute posterior predictive checks per group by using the T1 statistic (similar to chi square)
#' @param probit whether to use probit scale or probability scale
#' @seealso \code{\link{traitMPT}}
#' @author Daniel Heck
#' @export
getGroupMeans <- function(traitMPT, factor="all", probit=FALSE){

  if(is.null(traitMPT$mptInfo$predTable))
    stop("Model does not contain discrete predictors.")

  thetaNames <- traitMPT$mptInfo$thetaNames
  uniqueNames <- thetaNames[!duplicated(thetaNames[,"theta"]),"Parameter"]
  facLevelNames <- traitMPT$mptInfo$predFactorLevels
  if( all(factor == "all")){
    includeFactors <- names(facLevelNames)
  }else{
    if(any(! factor %in% names(facLevelNames)))
      stop("Factors defined by 'factor' are not factors in the fitted model!\n      Possible factors:  ",
           paste(names(facLevelNames)[!sapply(facLevelNames, is.null)], collapse=", "))
   includeFactors <- factor
    facLevelNames <- facLevelNames[factor]
  }
  S <- length(uniqueNames)

  summaryMat <- matrix(NA,0, 5)
  cnt <- 0
  for(s in 1:S){

    # parameter s of S
    muPosterior <- traitMPT$mcmc$BUGSoutput$sims.list[["mu"]][,s,drop=FALSE]

    # select only relevant predictors
    predTable <- traitMPT$mptInfo$predTable
    predTable <- subset(predTable,
                        predTable$prefix == "factor" &
                          predTable$theta == s &
                          predTable$Covariate %in% includeFactors)

    if(nrow(predTable) >0){
      factors <- unique(as.character(predTable$Covariate))
      combinations <-  expand.grid(facLevelNames[factors])
      combNames <- sapply(combinations, as.character)

      # get posterior samples for differences due to factors
      facPosterior <- list()
      for(j in 1:length(factors)){
        parnam <- as.character(subset(predTable, predTable$Covariate == factors[j])[,"covPar"])
        facPosterior[[j]] <- traitMPT$mcmc$BUGSoutput$sims.list[[parnam, drop=FALSE]]
        colnames(facPosterior[[j]]) <-  facLevelNames[[factors[j]]]
      }

      # splitMat <- unique(traitMPT$mptInfo$covData[,factors, drop=FALSE])
      for(i in 1:nrow(combinations)){
        samples <- muPosterior
        for(j in 1:ncol(combinations)){
          samples <- samples + facPosterior[[j]][,combinations[i,j]]
        }
        if(!probit)
          samples <- pnorm(samples)

        summaryMat <- rbind(summaryMat, c(Mean = mean(samples), SD = sd(samples),
                                          quantile(samples, c(.05,.5,.95))))
        rownames(summaryMat)[cnt <- cnt+1] <- paste0(uniqueNames[s],"_", paste0(paste0(
                                         colnames(combinations), "[",combNames[i,], "]"), collapse="_"))

      }
    }

  }
  summaryMat
}
