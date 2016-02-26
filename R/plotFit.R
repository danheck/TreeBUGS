
#' Plot Goodness of Fit
#'
#' Plots observed mean frequencies against sampled mean frequencies.
#'
#' @param fittedModel fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param ... arguments passed to \code{\link{boxplot}}
#' @export
plotFit <- function(fittedModel,...){

  dat <- fittedModel$mptInfo$dat
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  nam <- paste("response", TreeNames, "pred.mean", sep=".")
  select <- c(sapply(nam, grep, x=varnames(fittedModel$mcmc$mcmc)))
  # pred <- do.call("cbind", fittedModel$mcmc$BUGSoutput$sims.list[nam])
  pred <- do.call("rbind", fittedModel$mcmc$mcmc[,select])
  boxplot(pred, xaxt="n", main="Observed (red) and predicted (boxplot) mean frequencies", ...)
  axis(1, 1:ncol(dat), labels = colnames(dat))

  xx <- by(1:length(tree), tree, mean)
  axis(1, xx,  TreeNames, tick=F, line=NA, mgp=c(3, 2.5, 0))
  points(1:ncol(dat), colMeans(dat), col="red", cex=1.3, pch=19)
  abline(v=cumsum(table(tree))[1:(length(TreeNames)-1)]+.5)
}



