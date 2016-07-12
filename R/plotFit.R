
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

  # pred.mean <- colMeans(pred)
  # pred.95 <- apply(pred, 2, quantile, probs = c((1-ci)/2, 1- (1-ci)/2))
  # maxi <- max(pred.mean, pred.95, colMeans(dat))
  #
  # # Plot:
  # plot(pred.mean, pch=15, col="gray30",ylim=c(0, maxi),xaxt="n",cex=1, las=1,
  #      xlab="", ylab = "Mean frequencies",
  #      main="Observed (red) and predicted (gray) mean frequencies", ...)
  # segments(x0 = 1:ncol(pred.95), y0 = pred.95[1,], y1 = pred.95[2,], lwd = 2)
  # axis(1, 1:ncol(dat), labels = colnames(dat))
  # xx <- by(1:length(tree), tree, mean)
  # axis(1, xx,  TreeNames, tick=F, line=NA, mgp=c(3, 2.5, 0))
  # points(1:ncol(dat)+.1, colMeans(dat), col="red", cex=1, pch=16)
  # abline(v=cumsum(table(tree))[1:(length(TreeNames)-1)]+.5)
}



