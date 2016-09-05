#'	Read multiTree files
#'
#' Function to import MPT models from standard .eqn model files as used, for instance, by multiTree (Moshagen, 2010).
#'
#' @param file The (full path to the) file that specifies the multitree MPT file
#' @param restrictions Optional: The (full path to the) file that specifies which parameters should be constants and which should be equal. Alternatively: a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}
#' @param paramOrder if TRUE, the order of MPT parameters as interally used is printed.
#'
#' @details The file format should adhere to the standard .eqn-syntax (note that the first line is skipped and can be used for comments). In each line, a separate branch of the MPT model is specified using the tree label, category label, and the model equations in full form (multiplication sign `*` required; not abbreviations such as `a^2` allowed).
#'
#' As an example, the standard two-high threshold model (2HTM) is defined as follows:
#'
#'  \tabular{lllll}{
#' \code{Target } \tab \tab \code{Hit}             \tab \tab \code{Do} \cr
#' \code{Target}  \tab \tab \code{Hit}             \tab \tab \code{(1-Do)*g} \cr
#' \code{Target}  \tab \tab \code{Miss}            \tab \tab \code{(1-Do)*(1-g)} \cr
#' \code{Lure}    \tab \tab \code{FalseAlarm}      \tab \tab \code{(1-Dn)*g}  \cr
#' \code{Lure}    \tab \tab \code{CorrectReject}   \tab \tab \code{(1-Dn)*(1-g)} \cr
#' \code{Lure}    \tab \tab \code{CorrectReject  } \tab \tab \code{Dn}
#' }
#'
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")
#' readEQN(file = EQNfile, paramOrder = TRUE)
#'
#' # with equality constraint:
#' readEQN(file = EQNfile, restrictions = list("Dn=Do", "g=0.5"), paramOrder = TRUE)
#' @author Daniel Heck, Denis Arnold, Nina Arnold
#' @references Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. Behavior Research Methods, 42, 42-54.
#' @export
readEQN <- function(file, restrictions=NULL, paramOrder = FALSE){

  multiTreeDefinition = read.csv(file, header=F,
                                 blank.lines.skip = TRUE, sep= "",
                                 stringsAsFactors=F, skip = 1)	#read file

  # number of branches implied by number of rows in model file:
  numberOfBranches <- nrow(multiTreeDefinition)
  cols <- ncol(multiTreeDefinition)
  Tree <- data.frame(Tree = paste0("T_",multiTreeDefinition$V1),
                     Category = multiTreeDefinition$V2,
                     Equation = multiTreeDefinition$V3)
  Tree$Equation <- apply(multiTreeDefinition[,3:cols, drop=FALSE], 1, paste0, collapse="")

  TreeRestr <- thetaHandling(Tree, restrictions)

  allParameters <- getParameter(Tree)
  freeParameters <- subset(TreeRestr$SubPar, !duplicated(TreeRestr$SubPar$theta))$Parameter

  S <- length(freeParameters)
  numCat <- length(unique(Tree$Category))
  numTree <- length(unique(Tree$Tree))
  tt <- 1:S
  names(tt) <- freeParameters
  if(paramOrder){
    cat("Free parameters are used in the following order:\n")
    print(tt)
    cat("\n")

    if(!is.null(restrictions)){
      cat("Parameter constraints:\n")
      for(i in 1:length(TreeRestr$restrictions)){
        cat(TreeRestr$restrictions[[i]],"\n")
      }
      cat("\n")
    }

    if(S> numCat-numTree){
      cat("Note that the model is not identified and requires at least ",
          S-numCat+numTree ,"equality constraint(s).\n\n")
    }
  }

  # check MPT model
  par <- runif(length(allParameters))
  names(par) <- allParameters
  prob <- sapply(Tree$Equation, function(ff) eval(parse(text=ff), as.list(par)))
  sumPerTree <- as.vector(by(prob, Tree$Tree, sum))
  if(any(prob<0) | any(prob>1)){
    error <- paste0("Check .eqn-file. Model equations return values outside the iterval [0,1]:\n  ",
                    paste0("Line ", (1:length(prob))[prob<0 | prob>1],": ",
                           unique(Tree$Equation)[prob<0 | prob>1], collapse=", "))
    warning(error)
  }
  if(any(round(sumPerTree,8) != 1)){
    error <- paste0("Check .eqn-file. Probabilities do not sum up in trees:\n  ",
                    paste0(unique(Tree$Tree)[round(sumPerTree,8) != 1], collapse=", "))
    warning(error)
  }

  return(Tree)
}

# identifiability check
isIdentifiable <- function(S, Tree){
  numCat <- length(unique(Tree$Category))
  numTree <- length(unique(Tree$Tree))

  if(S> numCat-numTree){
    error <- paste0("Note that the model is not identified and requires at least ",
                    S-numCat+numTree ,"equality constraints.")
    warning(error)
  }
}
