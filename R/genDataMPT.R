#' Generate MPT Frequencies
#'
#' Uses a matrix of individual MPT parameters to generate MPT frequencies.
#' @param theta matrix of MPT parameters (rows: individuals; columns: parameters). Parameters are assigned by column names of the matrix. all of the parameters in the model file need to be included.
#' @param numItems number of responses per tree (a named vector with tree labels)
#' @inheritParams betaMPT
#' @param warning whether to show warning in case the naming of arguments does not match
#' @seealso \code{\link{genTraitMPT}} and \code{\link{genBetaMPT}} to generate data for latent normal/beta hierarchical distributions.
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")
#' theta <- matrix(c(.8,.4,.5,
#'                   .6,.3,.4), nrow=2, byrow=TRUE,
#'                 dimnames=list(NULL, c("Do","Dn","g")))
#' genDat <- genMPT(theta, c(Target=250, Lure=250),
#'                 EQNfile)
#' genDat
#' @export
genMPT <- function(theta, numItems, eqnfile, restrictions, warning=TRUE){
  if(missing(restrictions))
    restrictions <- NULL
  # read EQN
  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  Tree.restr <- thetaHandling(mergedTree,restrictions)
  # thetaNames <- getParameter(mergedTree)
  thetaNames <-  Tree.restr$SubPar[,1:2]
  thetaNames <- thetaNames[rownames(unique(thetaNames[2])),]$Parameter
  treeLabels <- unique(mergedTree$Tree)

  # get number of parmaeters/number of participants
  S <- length(thetaNames)
  if(is.vector(theta))
    theta <- matrix(theta, 1, dimnames=list(NULL, names(theta)))
  N <- nrow(theta)

  ################### check input + default values
  if(is.null(colnames(theta))){
    if(warning)
      warning("Colnames for theta are missing. Parameters are assigned by default as:\n  ",
              paste(thetaNames, collapse=", "))
    colnames(theta) <- thetaNames
  }
  if(is.null(names(numItems))){
    if(warning)
      warning("Tree labels for numitems are missing. Tree labels are assigned by default as:\n  ",
              paste(treeLabels, collapse=", "))
    names(numItems) <- treeLabels
  }else{
    names(numItems) <- paste0("T_", names(numItems))
  }
  theta <- checkThetaNames(theta, thetaNames)
  numItems <- checkNumItems(numItems, treeLabels)


  # matrix of response frequencies
  freq <- matrix(NA, N, nrow(mergedTree),
                 dimnames=list(NULL, mergedTree$Category))
  # consts <- Tree.restr$constants$sub
  # if(length(consts)>0)
  #   names(consts) <- Tree.restr$constants$Parameter
  colnames(theta) <- paste0("theta[",1:S,"]")
  eq <- gsub(",n","", Tree.restr$mergedTree$Equation, fixed=TRUE)
  for(n in 1:N){
    mergedTree$prob <- sapply(eq, function(ff) eval(parse(text=ff),
                                                    list(theta=theta[n,])))

    numTrees <- length(unique(mergedTree$Tree))
    for(k in 1:numTrees){
      sel <- mergedTree$Tree %in% names(numItems)[k]
      # cat <- findInterval(runif(numItems[k]), cumsum(mergedTree$prob[sel]))+1
      # catLabel <- mergedTree$Category[sel][cat]
      freq[n,mergedTree$Category[sel]] <- rmultinom(1, size = numItems[k],
                                                    prob=mergedTree$prob[sel])
      # table(factor(catLabel,
      #                                                levels=mergedTree$Category[sel]),
      #                                         exclude=NA)
    }
  }

  return(freq)
}
