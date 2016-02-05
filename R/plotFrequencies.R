

#' Plot Raw Frequencies
#'
#' Plots
#' @param fittedModel fitted hierarchical MPT model (see \code{\link{traitMPT}} or \code{\link{betaMPT}})
#' @param freq whether to plot absolute frequencies or relative frequencies (which sum up to one within each tree)
#' @export
plotFreq <- function(fittedModel, freq=TRUE){

  dat <- fittedModel$mptInfo$data

  # absolute frequencies
  if(freq){
    means <- colMeans(dat)
    dd <- dat[i,]
  }else{
    # relative frequencies (per tree)


  }


  K <- ncol(dat)
  plot(1:K, means, type="b", ylim=c(0, max(dat)), col=2, lwd=2)
  for(i in 1:K){
    lines(1:K, dd, col=adjustcolor(1, alpha=.6))
  }


}
