

#' Generate MPT Frequencies
#'
#' Uses a matrix of individual MPT parameters to generate MPT frequencies.
#'
#' @param theta matrix of MPT parameters (rows: individuals; columns: parameters). Parameters are assigned by column names of the matrix. all of the parameters in the model file need to be included.
#' @param numItems number of responses per tree (a named vector with tree labels)
#' @param eqnfile path to EQN file specifying the MPT model
#' @seealso \code{\link{genTraitMPT}} and \code{\link{genBetaMPT}} to generate data for latent normal/beta hierarchical distributions.
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")
#' theta <- matrix(c(.8,.4,.5, .6,.3,.4),
#'                 nrow=2, byrow=TRUE,
#'                 dimnames=list(NULL, c("Do","Dn","g")))
#' genDat <- genMPT(theta, c(Target=250, Lure=250),
#'                 EQNfile)
#' genDat
#' @export
genMPT <- function(theta, numItems, eqnfile){

  # read EQN
  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  thetaNames <- getParameter(mergedTree)
  treeLabels <- unique(mergedTree$Tree)

  # get number of parmaeters/number of participants
  S <- length(thetaNames)
  N <- nrow(theta)

  ################### check input + default values
  if(is.null(colnames(theta))){
    warning("Colnames for theta are missing. Parameters are assigned by default as:\n  ",
            paste(thetaNames, collapse=", "))
    colnames(theta) <- thetaNames
  }
  if(is.null(names(numItems))){
    warning("Tree labels for numitems are missing. Tree labels are assigned by default as:\n  ",
            paste(treeLabels, collapse=", "))
    names(numItems) <- treeLabels
  }else{
    names(numItems) <- paste0("T_", names(numItems))
  }
  checkThetaNames(theta, thetaNames)
  checkNumItems(numItems, treeLabels)


  # matrix of response frequencies
  freq <- matrix(NA, N, nrow(mergedTree),
                 dimnames=list(NULL, mergedTree$Category))
  for(n in 1:N){
    mergedTree$prob <- sapply(mergedTree$Equation,
                              function(ff) eval(parse(text=ff),
                                                as.list(theta[n,])))

    numTrees <- length(unique(mergedTree$Tree))
    for(k in 1:numTrees){
      sel <- mergedTree$Tree %in% names(numItems)[k]
      # cat <- findInterval(runif(numItems[k]), cumsum(mergedTree$prob[sel]))+1
      # catLabel <- mergedTree$Category[sel][cat]
      freq[n,mergedTree$Category[sel]] <- rmultinom(1, size = numItems[k], prob=mergedTree$prob[sel])
        # table(factor(catLabel,
        #                                                levels=mergedTree$Category[sel]),
        #                                         exclude=NA)
    }
  }

  return(freq)
}
