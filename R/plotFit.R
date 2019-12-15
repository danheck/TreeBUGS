#' Plot Posterior Predictive Mean Frequencies
#'
#' Plots observed means/covariances of individual frequencies against the means/covariances sampled from the posterior distribution (posterior predictive distribution).
#'
#' @inheritParams posteriorPredictive
#' @param stat whether to plot mean frequencies (\code{"mean"}) or covariances of individual frequencies (\code{"cov"})
#' @param ... arguments passed to \code{\link{boxplot}}
#'
#' @details If posterior predictive p-values were computed when fitting the model
#' (e.g., by adding the argument \code{traitMPT(...,ppp=1000)} ), the stored posterior samples are re-used for plotting.
#' Note that the last category in each MPT tree is dropped, because one category per multinomial distribution is fixed.
#'
#' @examples
#' \dontrun{
#' # add posterior predictive samples to fitted model:
#' fittedModel$postpred$freq.pred <-
#'      posteriorPredictive(fittedModel, M=1000)
#'
#' # plot model fit
#' plotFit(fittedModel, stat = "mean")
#' }
#' @export
plotFit <- function(fittedModel, M=1000, stat = "mean", ...){

  stat <- match.arg(stat, c("mean", "cov"))

  # get information about model:
  tree <- fittedModel$mptInfo$MPT$Tree
  cats <- fittedModel$mptInfo$MPT$Category
  dat <- fittedModel$mptInfo$dat[,cats]
  TreeNames <- unique(tree)

  # free categories (drop last category per tree):
  free_cats <- unlist(tapply(X = cats, INDEX = tree,
                             FUN = function(cat) cat[-length(cat)]))

  # get posterior predictive:
  if(is.null(fittedModel$postpred) | M != 1000){
    freq.list <- posteriorPredictive(fittedModel, M=M)
  }else{
    freq.list <- fittedModel$postpred$freq.pred
  }

  if(stat == "mean"){
    # Plot mean frequencies:

    pred <- t(sapply(freq.list, colMeans))
    boxplot(pred[,free_cats], xaxt="n", col="gray",
            main="Observed (red) and predicted (boxplot) mean frequencies", las=1, ...)
    axis(1, 1:length(free_cats), labels = free_cats)
    xx <- by(seq_along(free_cats), names(free_cats), mean)
    axis(1, xx,  TreeNames, tick=F, line=NA, mgp=c(3, 2.5, 0))
    points(1:length(free_cats), colMeans(dat)[free_cats],
           col="red", cex=1.4, pch=17)
    abline(v = cumsum(table(tree) - 1)[1:(length(TreeNames)-1)]+.5, col="gray")

  } else if (stat == "cov"){
    # Plot covariance of frequencies:

    nams <- outer(free_cats, free_cats, paste, sep="-")
    sel_cov <- nams[upper.tri(nams, diag = TRUE)]
    K <- length(sel_cov)

    # observed/predicted
    c.obs <- cov(dat[,free_cats])
    c.pred <- sapply(freq.list, function(xx){
      cc <- cov(xx[,free_cats])
      cc[upper.tri(cc, diag=TRUE)]
    })

    boxplot(t(c.pred), col="gray", ylab="Covariance",
            main="Observed (red) and predicted (gray) covariances",
            xaxt="n", las=1, ...)
    abline(h=0, lty=1, col="gray")
    axis(1, 1:K, labels = nams[upper.tri(nams, diag=TRUE)], las=2)
    points(1:K, c.obs[upper.tri(c.obs, diag=TRUE)], col=2, pch=17)
    abline(v = cumsum(seq(nrow(c.obs), 2, -1))+.5, col="lightgray")
  }
}



