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
#' EQNfile <- paste0(.libPaths()[1], "/TreeBUGS/MPTmodels/2htm.eqn")
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
  Tree <- data.frame(Tree = multiTreeDefinition$V1,
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

# #' Unique branches for each Tree with the same answer and sum the corresponing formulas
# #'
# #' @param Tree Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
mergeBranches<-function(Tree){ # OLD ,DataNames){ # Unique branches for each Tree with one answer and sum the corresponing formulas

  treeNames <- sort(unique(Tree$Tree))
  catNames <- as.list(by(Tree$Category, Tree$Tree, function(xx) sort(unique(xx))))
  names(catNames) <- treeNames

#   for(i in 1:length(Tree)){
#     TreesUniqueAnswers[[i]]=unique(Tree[Tree$Tree==Tree[i],]$Category)
#   }

  NewTree <- data.frame(Tree = rep(treeNames, sapply(catNames, length)),
                        Category = unlist(catNames),
                        Equation = "...", stringsAsFactors = FALSE)
  for(tt in 1:length(treeNames)){
    for(cc in 1:length(catNames[[tt]])){
      selTree <- Tree$Tree == treeNames[tt] & Tree$Category == catNames[[tt]][cc]
      selNewTree <- NewTree$Tree == treeNames[tt] & NewTree$Category == catNames[[tt]][cc]
      NewTree$Equation[selNewTree] <- paste(Tree[selTree,]$Equation, collapse="+")
    }
  }

#   NewTree <- data.frame(Tree="character",
#                         Category="character",
#                         Equation="character",
#                         stringsAsFactors=F)
#   NewTreesFormula=TreesUniqueAnswers
#
#   for(i in 1:length(Tree)){
#     for(j in 1:length(TreesUniqueAnswers[[i]])){
#       NewTreesFormula[[i]][j] <- paste(TreeData$Equation[TreeData$Tree==Tree[i] &
#                                                            TreeData$Category==TreesUniqueAnswers[[i]][j]],collapse="+")
#     }
#   }

#   for(i in 1:length(Tree)){
#     for(j in 1:length(TreesUniqueAnswers[[i]])){
#       NewTrees=rbind(NewTrees,data.frame(Tree=Tree[i],
#                                          Category=TreesUniqueAnswers[[i]][j],
#                                          Equation=NewTreesFormula[[i]][j],stringsAsFactors=F))
#     }
#   }
#
#   mergedTree=NewTrees[-1,]
#   mergedTree=mergedTree[order(mergedTree$Tree),]

  return(NewTree)
}



# #' Extract all parameter from the formulas of a given model
# #'
# #' @param TreeData Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
getParameter<-function(TreeData){

  Parameter=unique(unlist(strsplit(TreeData$Equation,
                                   split="\\*|\\(|\\)|\\-|\\+")))
  r=c(which(nchar(Parameter)==0), grep("^[0-9]+$|^[0-9]+\\.[0-9]+", Parameter))
  Parameter=Parameter[-r]
  Parameter=c(sort(Parameter[grep("[A-Z]",Parameter)]),
              sort(Parameter[-grep("[A-Z]",Parameter)]))
  return(sort(Parameter))

}


# #' Read the subject data from file
# #'
# #' @param data the data as data.frame
# #' @param Category is the unique of the $Category from TreeData which is $Category after mergeBrachnes()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
readSubjectData<-function(data,Category){

  if(sum(names(data) %in% Category)!=length(Category)){
    if(dim(data)[2]!=length(Category)){
      cat("Number of categories differs in eqn and csv file")
      return(-1)
    }
    else{
      cat("At least one name of the categories differs in eqn and csv file")
      return(-1)
    }
  }

  return(data)
}


# #'  Set constants, replace parameters for thetas
# #'
# #' @param mergedTree is the output of mergeBranches
# #' @param restrictions Either the (full path to the) file that specifies which parameters should be constants and which should # be equal; or a list of parameter restrictions
# #' @author Nina R. Arnold, Denis Arnold, Daniel Heck
# #' @export
thetaHandling <-function(mergedTree, restrictions){

  Parameter=getParameter(mergedTree)

  SubPar <- data.frame(Parameter = Parameter,
                       theta = 1:length(Parameter),
                       sub="", stringsAsFactors=F)

  ############ only if restrictions are included:
  if(!is.null(restrictions)){

    # restrictions given as a list
    if(is.list(restrictions)){
      restrVector <- as.vector(unlist(restrictions))

    # restrictions given as a model file
    }else{
      restrVector <- read.csv(restrictions, header=F,stringsAsFactors=F)$V1
      restrictions <- as.list(restrVector)
    }
    restrVector <- gsub(" ", "", restrVector)

    for(k in 1:length(restrVector)){
      splitRestr <- strsplit(restrVector[k], "=")[[1]]
      if(length(splitRestr) == 1){
        warning("Restriction not well defined: Equality sign '=' missing in:\n  ",splitRestr)
      }else{
        index <- match(splitRestr, SubPar$Parameter)
        suppressWarnings(consts <- as.numeric(splitRestr))
        # consts <- consts[!is.na(consts)]

        if(sum(!is.na(consts)) == 0){
          # only parameters without constants
          if(any(is.na(index))){
            error <- paste0("Restriction contains parameters not contained in the model:\n  ",
                            paste(splitRestr, collapse="="))
            stop(error)
          }
          # replace index
          SubPar[index[2:length(index)], "theta"]<- index[1]
        }else if(sum(!is.na(consts)) == 1){
          # contrained to constant values
          CONST <- consts[!is.na(consts)]
          if(CONST <0 | CONST >1){
            error <- paste0("Check parameter restrictions. Constants are not in the interval [0,1]: ",
                            restrVector[k])
            warning(error)
          }
          SubPar[index[!is.na(index)], "theta"] <- - CONST
        }else{
          stop("Restrictions should not contain more than one constant!")
        }
      }
    }
  }

  isConstant <- SubPar$theta <= 0

  ############## replaced constant parameter values:
  if(sum(isConstant) > 0){
    constants <- SubPar[isConstant,, drop=FALSE]
    SubPar <- SubPar[!isConstant,, drop=FALSE]
    constants$sub <- - constants$theta

    for(j in 1:nrow(constants)){
      mergedTree$Equation <- gsub(
        pattern = paste0("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
                         constants$Parameter[j],
                         "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
        replacement = paste0("\\1",constants$sub[j],"\\2"),
        mergedTree$Equation, perl=T)
    }
#     }
  }else{
    constants <- NULL
  }

  # use new, increasing indices 1....S for parameters:
  tmp <- unique(SubPar$theta)
  for(tt in 1:length(tmp)){
    SubPar$theta[SubPar$theta == tmp[tt]] <- tt
  }

  ############## replaced constant parameter values:
  SubPar$sub <- paste0("theta[", SubPar$theta, ",n]")
  for(j in 1:nrow(SubPar)){
    mergedTree$Equation <- gsub(
      pattern = paste0("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
                       SubPar$Parameter[j],
                       "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
      replacement = paste0("\\1",SubPar$sub[j],"\\2"),
      mergedTree$Equation, perl=T)
  }

  #################### OLD VERSION: #############################

#   for(i in 1:nrow(SubPar)){
#     mergedTree$Equation <- gsub(pattern = SubPar$Parameter[i],
#                                 replacement = SubPar$sub[i],
#                                 mergedTree$Equation, fixed = TRUE)
#   }

#   partrack=Parameter
#
#   # Handle the constants
#   if(missing(restrictions) | is.null(restrictions)){
#     # no constraints: TODO
#     SubPar <- NULL
#     numberSubs <- 0
#   }else{
#     # with constraints:
#
#     if(is.list(restrictions)){
#       # restrictions given as a list
#       SubPar <- data.frame(V1 = as.vector(unlist(restrictions)),
#                            stringsAsFactors = FALSE)
#     }else{
#       # restrictions given as a model file
#       SubPar=read.csv(restrictions,header=F,stringsAsFactors=F)
#     }
#
#     # actual replacement of parameters:
#     numberSubs <- grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1) # grep("=", SubPar$V1) #
#   }
#   replaceConst=data.frame(const="character",sub="character",stringsAsFactors=F)
#   if(any(numberSubs > 0)){
#     for(i in 1:length(numberSubs)){
#       X=unlist(strsplit(SubPar$V1[numberSubs[i]],"="))
#       number=grep("^[0-9]+$|^[0-9]+\\.[0-9]+$",X)
#       if(length(number)>1){
#         cat("Variable set to multiple constants")
#         return(-1)
#       }
#       else{
#         for(i in 1:length(X)){
#           if(i!=number){
#             replaceConst=rbind(replaceConst,
#                                data.frame(const=X[i], sub=X[number]))
#           }
#         }
#
#       }
#     }
#     replaceConst = replaceConst[-1,]
#     for(i in 1:dim(replaceConst)[1]){
#       mergedTree$Equation=gsub(paste0(
#         "(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
#         replaceConst$const[i],
#         "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
#         paste0("\\1", replaceConst$sub[i],"\\2"),
#         mergedTree$Equation,perl=T)
#     }
#   }
#
#   partrack=partrack[-which(partrack%in%replaceConst$const)]
#
#   # Handle parameters set equal
#   if(length(grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1))>0){
#     SubPar=SubPar[-grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1),]
#   }
#   SubPar=strsplit(SubPar,split="=")
#
#   theta=numeric()
#
#   for(i in 1:length(SubPar)){
#     for(j in 1:length(SubPar[[i]])){
#       theta=c(theta,i)
#     }
#
#   }
#   SubPar=data.frame(Par=unlist(SubPar),
#                     theta=theta,
#                     stringsAsFactors=F)
#   SubPar$sub=paste("theta[",SubPar$theta,",n]",sep="")
#
#   partrack=partrack[-which(partrack%in%SubPar$Par)]
#   if(length(partrack)>0){
#     for(i in 1:length(partrack)){
#       SubPar=rbind(SubPar,data.frame(Par=partrack[i],
#                                      theta=max(SubPar$theta+1),
#                                      sub=paste("theta[",max(SubPar$theta)+1,",n]",sep="")))
#     }
#   }
#
  # for(j in 1:length(SubPar$Par)){


  output=list(SubPar = SubPar, mergedTree = mergedTree,
              constants=constants, restrictions = restrictions)

  return(output)

}
