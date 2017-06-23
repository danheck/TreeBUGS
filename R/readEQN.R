#'	Read multiTree files
#'
#' Function to import MPT models from standard .eqn model files as used, for instance, by multiTree (Moshagen, 2010).
#'
#' @param file The (full path to the) file that specifies the MPT model (standard .eqn syntax). Note that category labels must start with a letter (different to multiTree) and match the column names of \code{data}. Alternatively, the EQN-equations can be provided within R as a character value (see examples). Note that the first line of an .eqn-file is reserved for comments and always ignored.
#' @inheritParams betaMPT
#' @param paramOrder if TRUE, the order of MPT parameters as interally used is printed.
#' @param parse whether to return a parsed MPT model description in terms of the matrices \eqn{a} and \eqn{b} (the powers of the \eqn{\theta} and \eqn{(1-\theta)}, respectively, and the vector of constants \eqn{c}. Each branch probability is then given as \eqn{c_{i}  \prod_{s} \theta^{a_{i,s}}(1-\theta)^{b_{i,s}})}
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
#' EQNfile <- system.file("MPTmodels/2htm.eqn",
#'                        package="TreeBUGS")
#' readEQN(file = EQNfile, paramOrder = TRUE)
#'
#' # with equality constraint:
#' readEQN(file = EQNfile, restrictions = list("Dn=Do", "g=0.5"),
#'         paramOrder = TRUE)
#'
#' # define MPT model directly within R
#' model <-
#'   "2-High Threshold Model (2HTM)
#'   old hit d
#'   old hit (1-d)*g
#'   old miss (1-d)*(1-g)
#'   new fa (1-d)*g
#'   new cr (1-d)*(1-g)
#'   new cr d"
#' readEQN(model, paramOrder=TRUE)
#' @author Daniel Heck, Denis Arnold, Nina Arnold
#' @references Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. Behavior Research Methods, 42, 42-54.
#' @export
readEQN <- function(file,
                    restrictions=NULL,
                    paramOrder = FALSE,
                    parse=FALSE){
  if(missing(restrictions)) restrictions <- NULL

  isPath <- !grepl("\n", x=file, ignore.case = TRUE)
  if(!isPath){
    model <- file
    file <- tempfile(pattern = "MPTmodel", tmpdir = tempdir(), fileext = ".eqn")
    cat(paste0(model,"\n"), file=file)
  }
  # read first line if it contains model equations
  multiTreeDefinition <- read.csv(file, header=F, comment.char = "#",
                                  blank.lines.skip = TRUE, sep= "",
                                  stringsAsFactors=F, skip = 1)


  # number of branches implied by number of rows in model file:
  numberOfBranches <- nrow(multiTreeDefinition)
  cols <- ncol(multiTreeDefinition)
  Tree <- data.frame(Tree = paste0("T_",multiTreeDefinition$V1),
                     Category = multiTreeDefinition$V2,
                     Equation = multiTreeDefinition$V3)
  Tree$Equation <- apply(multiTreeDefinition[,3:cols, drop=FALSE], 1, paste0, collapse="")

  TreeRestr <- thetaHandling(Tree, restrictions)

  allParameters <- getParameter(Tree)
  suppressWarnings(
    free <- !duplicated(TreeRestr$SubPar$theta) &
      is.na(as.numeric(TreeRestr$SubPar$Parameter))
  )
  freeParameters <- subset(TreeRestr$SubPar, free)$Parameter

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
                    paste0(unique(Tree$Tree)[round(sumPerTree,8) != 1], collapse=", "),
                    "\n  (note that the first line of the .eqn-file ist ignored!)")
    warning(error)
  }

  if(parse){
    tmp <- Tree
    mpt <- parseEQN(Tree)
    Tree <- parseRestrictions(mpt, restrictions)
    Tree$Table <- tmp
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


# mergeTree <- function(eqnfile, data, restrictions){
#
#   # MPT structure for JAGS
#   Tree <- readEQN(eqnfile)
#   mergedTree <- mergeBranches(Tree)
#
#   tHoutput <- thetaHandling(mergedTree,restrictions)
#   SubPar <- tHoutput$SubPar
#   mergedTree <- tHoutput$mergedTree
#   fixedPar <- tHoutput$fixedPar
# }
