
#' Plot Posterior Predictive Mean Frequencies
#'
#' Plots observed mean frequencies and boxplot of mean frequencies sampled from posterior distribution (posterior predictive).
#'
#' @inheritParams posteriorPredictive
#' @param ... arguments passed to \code{\link{boxplot}}
#'
#' @details If posterior predictive p-values were computed when fitting the model (e.g., \code{traitMPT(...,ppp=1000)} ), the stored posterior samples are re-used for plotting.
#' @export
plotFit <- function(fittedModel, M=1000, ...){

  # get information about model:
  dat <- fittedModel$mptInfo$dat
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)

  # get posterior predictive:
  if(is.null(fittedModel$postpred) | M != 1000){
    freq.list <- posteriorPredictive(fittedModel, M=M)
  }else{
    freq.list <- fittedModel$postpred$freq.pred
  }
  pred <- t(sapply(freq.list, colMeans))

  # Plot:
  boxplot(pred, xaxt="n",
          main="Observed (red) and predicted (boxplot) mean frequencies", ...)
  axis(1, 1:ncol(dat), labels = colnames(dat))
  xx <- by(1:length(tree), tree, mean)
  axis(1, xx,  TreeNames, tick=F, line=NA, mgp=c(3, 2.5, 0))
  points(1:ncol(dat), colMeans(dat), col="red", cex=1.4, pch=17)
  abline(v=cumsum(table(tree))[1:(length(TreeNames)-1)]+.5)
}



