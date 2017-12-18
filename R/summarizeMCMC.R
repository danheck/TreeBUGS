
#' MCMC Summary
#'
#' TreeBUGS-specific summary for \code{mcmc.list}-objects.
#' @param mcmc a \code{\link[coda]{mcmc.list}} object
#' @param batchSize size of batches of parameters used to reduce memory load when computing (univariate) Rhat statistics
#' @param probs quantile probabilities used to compute credibility intervals
#' @export
summarizeMCMC <- function(mcmc, batchSize=100, probs = c(.025,.50,.975)){
  if(class(mcmc) %in% c("traitMPT", "betaMPT","simpleMPT"))
    mcmc <- mcmc$runjags$mcmc
  if(class(mcmc) == "runjags")
    mcmc <- mcmc$mcmc

  mcmc.mat <- do.call("rbind", mcmc)
  summTab <- cbind("Mean"=apply(mcmc.mat,2,mean, na.rm = TRUE),
                   "SD"=apply(mcmc.mat,2,sd, na.rm = TRUE),
                   t(apply(mcmc.mat, 2, quantile, probs, na.rm = TRUE)),
                   "Time-series SE"=NA, "n.eff" = NA ,
                   "Rhat" = NA, "R_95%"=NA)
  #     summ[[1]][,1:2], summ[[2]], "Time-series SE"=summ[[1]][,4]
  rm(mcmc.mat)
  gc(verbose=FALSE)
  rn <- rownames(summTab)
  # sel.notT1 <- setdiff(1:nrow(summTab), union(grep("T1", rn), grep(".pred.mean", rn)))
  try({
    summTab[,"n.eff"] <- round(effectiveSize(mcmc))
    summTab[,"Time-series SE"] <- summTab[,"SD"] / sqrt(summTab[,"n.eff"]  )
  }, silent = TRUE)
  gc(verbose=FALSE)
  if(is.list(mcmc) && length(mcmc) > 1){
    try({
      # batchSize <- 200
      # split for 100 variables per batch
      n.summ <- nrow(summTab)
      if(n.summ >batchSize){
        for(ii in 1:(n.summ %/% batchSize)){
          idx <- (ii-1)*batchSize + 1:batchSize
          # summ.idx <- sel.notT1[idx]
          summTab[idx,c("Rhat", "R_95%")] <- gelman.diag(mcmc[,idx], multivariate=FALSE)[[1]]
          gc(verbose=FALSE)
        }
        if((ii*batchSize+1) < n.summ){
          # summ.idx <- sel.notT1[(ii*batchSize+1):n.summ]
          idx <- (ii*batchSize+1):n.summ
          summTab[idx,c("Rhat", "R_95%")] <- gelman.diag(mcmc[,idx], multivariate=FALSE)[[1]]
        }
      }else{
        summTab[,c("Rhat", "R_95%")] <- gelman.diag(mcmc, multivariate=FALSE)[[1]]
      }
    })}
  #   if(any(is.na(summTab[,"Rhat"])))
  #     warning("Gelman-Rubin convergence diagnostic Rhat could not be computed.")
  # n.eff <-
  gc(verbose=FALSE)

  summTab
}
