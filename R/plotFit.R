
#' Plot Posterior Predictive Mean Frequencies
#'
#' Plots observed mean frequencies and boxplot of mean frequencies sampled from posterior distribution (posterior predictive).
#'
#' @inheritParams posteriorPredictive
#' @param ... arguments passed to \code{\link{boxplot}}
#' @export
plotFit <- function(fittedModel, M=1000, ...){

  # get information about model:
  dat <- fittedModel$mptInfo$dat
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)

  # get posterior predictive:
  freq.list <- posteriorPredictive(fittedModel, M=M)
  pred <- t(sapply(freq.list, colMeans))

  # Plot:
  boxplot(pred, xaxt="n", main="Observed (red) and predicted (boxplot) mean frequencies", ...)
  axis(1, 1:ncol(dat), labels = colnames(dat))
  xx <- by(1:length(tree), tree, mean)
  axis(1, xx,  TreeNames, tick=F, line=NA, mgp=c(3, 2.5, 0))
  points(1:ncol(dat), colMeans(dat), col="red", cex=1.3, pch=19)
  abline(v=cumsum(table(tree))[1:(length(TreeNames)-1)]+.5)
}



