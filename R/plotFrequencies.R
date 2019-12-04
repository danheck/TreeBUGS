


#' Plot Raw Frequencies
#'
#' Plot observed individual and mean frequencies.
#'
#' @param x either a fitted hierarchical MPT model (see \code{\link{traitMPT}}, \code{\link{betaMPT}});
#'     or a matrix/data frame of response frequencies (can be provided as a path to a .csv-file with individual frequencies).
#' @param freq whether to plot absolute frequencies or relative frequencies
#'     (which sum up to one within each tree; only if \code{x} is a hierarchical model or if \code{eqnfile} is provided)
#' @param select a numeric vector with participant indices to select which raw frequencies to plot
#'     (default: \code{"all"})
#' @param boxplot if \code{FALSE}, lines and points are drawn instead of boxplots
#' @param eqnfile optional: EQN description of an MPT model, that is, either the path to an EQN file or as a character string
#'     (only used if \code{x} refers to a matrix/data frame or .csv-file)
#' @param ... further arguments passed to \code{boxplot} and \code{plot}
#' @export
#' @examples
#' # get frequency data and EQN file
#' freq <- subset(arnold2013, group == "encoding", select = -(1:4))
#' eqn <- system.file("MPTmodels/2htsm.eqn", package="TreeBUGS")
#' plotFreq(freq, eqnfile = eqn)
#' plotFreq(freq, freq = FALSE, eqnfile = eqn)
plotFreq <- function(x, freq=TRUE, select="all", boxplot=TRUE, eqnfile,...){


  if(inherits(x, c("betaMPT", "traitMPT"))){
    dat <- x$mptInfo$data
  }else if(inherits(x, "character")){
    dat <- read.csv(x)
  }else{
    try(dat <- as.data.frame(x))
  }

  if(inherits(x, c("betaMPT", "traitMPT"))){
    treeNames <- x$mptInfo$MPT$Tree
    treeLabels <- unique(treeNames)
  }else if(!missing(eqnfile)){
    tmp <- unique(readEQN(eqnfile)[,1:2])
    treeNames <- tmp$Tree
    treeLabels <- unique(treeNames)
    try(dat <- dat[,colnames(dat) %in% tmp$Category])
  }else{
    treeNames <- rep("", ncol(dat))
    treeLabels <- ""
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

  # absolute frequencies
  if(!freq){
    # relative frequencies (per tree)
    for(t in 1:length(treeLabels)){
      sel <- treeNames == treeLabels[t]
      dat[,sel] <- dat[,sel] / rep(rowSums(dat[,sel]), each=sum(sel))
    }
  }

  means <- colMeans(dat)

  if(boxplot == TRUE){
    boxplot(dat, ylab=ifelse(freq, "Absolute frequency", "Relative frequency (per tree)"),
            xlab="", main=ifelse(freq, "Absolute frequency", "Relative frequency (per tree)"), las=1, ...)
    lines(1:K, means, col="red", lwd=2)
  }else{
    plot(1:K, rep(NA, K), ylim=c(0, max(dat)), col=1, lwd=3, pch=16, xaxt="n", las=1, ...,
         ylab=ifelse(freq, "Absolute frequency", "Relative frequency (per tree)"),
         xlab="", main=ifelse(freq, "Absolute frequency", "Relative frequency (per tree)"))
    axis(1, 1:K, colnames(dat))
    for(treelab in treeLabels){
      sel <- treeNames == treelab
      for(i in 1:N){
        lines((1:K)[sel], dat[i,sel], col=rainbow(N, alpha=.4)[i])
        points((1:K)[sel], dat[i,sel], col=rainbow(N, alpha=.6)[i], pch=16)
      }
      lines((1:K)[sel], means[sel], col=1, lwd=3)
    }

  }

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
