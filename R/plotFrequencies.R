


#' Plot Raw Frequencies
#'
#' Plot observed raw frequencies
#' @param fittedModel fitted hierarchical MPT model (see \code{\link{traitMPT}} or \code{\link{betaMPT}}). Can also be a path to a .csv-file with individual frequencies or a matrix/data frame.
#' @param freq whether to plot absolute frequencies or relative frequencies (which sum up to one within each tree; only if \code{fittedModel} is a hierarchical model)
#' @param select a numeric vector with participant indices to select which raw frequencies to plot (default: \code{"all"})
#' @param eqnfile optional argument to get MPT tree structure if \code{fittedModel} is the path to a .csv-file or a matrix/data frame
#' @export
plotFreq <- function(fittedModel, freq=TRUE, select="all", eqnfile){


  if(class(fittedModel) %in% c("betaMPT", "traitMPT")){
    dat <- fittedModel$mptInfo$data
  }else if(class(fittedModel) == "character"){
    dat <- read.csv(fittedModel)
  }else{
    try(dat <- as.data.frame(fittedModel))
  }

  K <- ncol(dat)
  N <- nrow(dat)

  if(select == "all"){
    select <- 1:N
  }else{
    if(!is.numeric(select)  || any(select != round(select)) )
      stop("Please use an integer vector to select participants.")
    dat <- dat[select, ,drop=FALSE]
    N <- nrow(dat)
  }

  if(class(fittedModel) %in% c("betaMPT", "traitMPT")){
    treeNames <- fittedModel$mptInfo$MPT$Tree
    treeLabels <- unique(treeNames)
  }else if(!missing(eqnfile)){
    tmp <- unique(readEQN(eqnfile)[,1:2])
    treeNames <- tmp$Tree
    treeLabels <- unique(treeNames)
  }else{
    treeNames <- rep("", ncol(dat))
    treeLabels <- ""
  }

  # absolute frequencies
  if(!freq){
    # relative frequencies (per tree)
    for(t in 1:length(treeLabels)){
      sel <- treeNames == treeLabels[t]
      dat[,sel] <- dat[,sel] / rep(rowSums(dat[,sel]), each=sum(sel))
    }
  }

  means <- colMeans(dat)


  plot(1:K, means, ylim=c(0, max(dat)), col=1, lwd=3, pch=16, xaxt="n",
       ylab=ifelse(freq, "Absolute frequency", "Relative frequency (per tree)"),
       xlab="",
       main="Raw Frequencies")
  axis(1, 1:K, colnames(dat))
  for(i in 1:N){
    lines(1:K, dat[i,], col=rainbow(N, alpha=.5)[i])
    points(1:K, dat[i,], col=rainbow(N, alpha=.6)[i], pch=16)
  }
  lines(1:K, means, col=1, lwd=3)

  xt <- .5
  for(k in 2:K){
    if( treeNames[k] != treeNames[k-1]){
      abline(v=k-.5)
      xt <- c(xt, k-.5)
    }
  }
  xt <- c(xt, K+.5)
  axis(1, xt[1:(length(xt)-1)]+ diff(xt)/2,
       treeLabels, mgp=c(100,3,10))
}
