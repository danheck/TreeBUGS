#'	Read multiTree files
#'
#' Function to import MPT models from standard .eqn model files as used, for instance, by multiTree (Moshagen, 2010).
#'
#' @param file The (full path to the) file that specifies the multitree MPT file
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
#' @author Daniel Heck, Denis Arnold, Nina Arnold
#' @references Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. Behavior Research Methods, 42, 42-54.
#' @export
readEQN <- function(file, paramOrder = FALSE){

  multiTreeDefinition = read.csv(file, header=F,
                                 blank.lines.skip = TRUE, sep= "",
                                 stringsAsFactors=F, skip = 1)	#read file

  # number of branches implied by number of rows in model file:
  numberOfBranches <- nrow(multiTreeDefinition)
  cols <- ncol(multiTreeDefinition)
  TreeData <- data.frame(Tree = multiTreeDefinition$V1,
                         Category = multiTreeDefinition$V2,
                         Equation = multiTreeDefinition$V3)
  TreeData$Equation <- apply(multiTreeDefinition[,3:cols, drop=FALSE], 1, paste0, collapse="")

  nn <- getParameter(TreeData)
  S <- length(nn)
  numCat <- length(unique(TreeData$Category))
  numTree <- length(unique(TreeData$Tree))
  tt <- 1:S
  names(tt) <- nn
  if(paramOrder){
    cat("Parameters are used in the following order:\n")
    print(tt)
    cat("\n")

    if(S> numCat-numTree){
      cat("Note that the model is not identified and requires at least ",
          S-numCat+numTree ,"equality constraints.\n\n")
    }
  }

  # check MPT model
  par <- runif(S)
  names(par) <- nn
  prob <- sapply(TreeData$Equation, function(ff) eval(parse(text=ff), as.list(par)))
  sumPerTree <- as.vector(by(prob, TreeData$Tree, sum))
  if(any(prob<0) | any(prob>1)){
    error <- paste0("Check .eqn-file. Model equations return values outside the iterval [0,1]:\n  ",
                    paste0("Line ", (1:length(prob))[prob<0 | prob>1],": ",
                           unique(TreeData$Equation)[prob<0 | prob>1], collapse=", "))
    warning(error)
  }
  if(any(round(sumPerTree,8) != 1)){
    error <- paste0("Check .eqn-file. Probabilities do not sum up in trees:\n  ",
                    paste0(unique(TreeData$Tree)[round(sumPerTree,8) != 1], collapse=", "))
    warning(error)
  }

  return(TreeData)
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
# #' @param TreeData Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
mergeBranches<-function(Tree){ # OLD ,DataNames){ # Unique branches for each Tree with one answer and sum the corresponing formulas

  treeNames <- sort(unique(Tree$Tree))
  catNames <- as.list(by(Tree$Category, Tree$Tree, function(xx) sort(unique(xx))))
  names(catNames) <- treeNames

#   for(i in 1:length(Tree)){
#     TreesUniqueAnswers[[i]]=unique(TreeData[TreeData$Tree==Tree[i],]$Category)
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
#   MergedTree=NewTrees[-1,]
#   MergedTree=MergedTree[order(MergedTree$Tree),]

  return(NewTree)
}



# #' Extract all parameter from the formulas of a given model
# #'
# #' @param TreeData Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
getParameter<-function(TreeData){

  Parameter=unique(unlist(strsplit(TreeData$Equation,split="\\*|\\(|\\)|\\-|\\+")))
  r=c(which(nchar(Parameter)==0),grep("^[0-9]+$|^[0-9]+\\.[0-9]+",Parameter))
  Parameter=Parameter[-r]
  Parameter=c(sort(Parameter[grep("[A-Z]",Parameter)]),sort(Parameter[-grep("[A-Z]",Parameter)]))
  return(Parameter)

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
# #' @param MergedTree is the output of mergeBranches
# #' @param restrictions Either the (full path to the) file that specifies which parameters should be constants and which should # be equal; or a list of parameter restrictions
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
thetaHandling<-function(MergedTree,restrictions){

  Parameter=getParameter(MergedTree)
  partrack=Parameter

  # Handle the constants
  if(missing(restrictions) | is.null(restrictions)){
    # no constraints: TODO
    SubPar=NULL

  }else{
    # with constraints:

    if(is.list(restrictions)){
      # restrictions given as a list
      SubPar <- data.frame(V1 = as.vector(unlist(restrictions)),
                           stringsAsFactors = FALSE)
    }else{
      # restrictions given as a model file
      SubPar=read.csv(restrictions,header=F,stringsAsFactors=F)
    }

    # actual replacement of parameters:

    numberSubs <- grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1) # grep("=", SubPar$V1) #
    replaceConst=data.frame(const="character",sub="character",stringsAsFactors=F)
    if(length(numberSubs)>0){
      for(i in 1:length(numberSubs)){
        X=unlist(strsplit(SubPar$V1[numberSubs[i]],"="))
        number=grep("^[0-9]+$|^[0-9]+\\.[0-9]+$",X)
        if(length(number)>1){
          cat("Variable set to multiple constants")
          return(-1)
        }
        else{
          for(i in 1:length(X)){
            if(i!=number){
              replaceConst=rbind(replaceConst,
                                 data.frame(const=X[i], sub=X[number]))
            }
          }

        }
      }
      replaceConst = replaceConst[-1,]
      for(i in 1:dim(replaceConst)[1]){
        MergedTree$Equation=gsub(paste0(
          "(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
          replaceConst$const[i],
          "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
          paste0("\\1", replaceConst$sub[i],"\\2"),
          MergedTree$Equation,perl=T)
      }
    }

    partrack=partrack[-which(partrack%in%replaceConst$const)]

    # Handle parameters set equal
    if(length(grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1))>0){
      SubPar=SubPar[-grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1),]
    }
    SubPar=strsplit(SubPar,split="=")

    theta=numeric()

    for(i in 1:length(SubPar)){
      for(j in 1:length(SubPar[[i]])){
        theta=c(theta,i)
      }

    }
    SubPar=data.frame(Par=unlist(SubPar),
                      theta=theta,
                      stringsAsFactors=F)
    SubPar$sub=paste("theta[",SubPar$theta,",n]",sep="")

    partrack=partrack[-which(partrack%in%SubPar$Par)]
    if(length(partrack)>0){
      for(i in 1:length(partrack)){
        SubPar=rbind(SubPar,data.frame(Par=partrack[i],
                                       theta=max(SubPar$theta+1),
                                       sub=paste("theta[",max(SubPar$theta)+1,",n]",sep="")))
      }
    }

    for(j in 1:length(SubPar$Par)){
      MergedTree$Equation=gsub(paste("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",SubPar$Par[j],"(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)",sep=""),
                               paste("\\1",SubPar$sub[j],"\\2",sep=""),
                               MergedTree$Equation,perl=T)
    }
  }

  output=list(SubPar,MergedTree)

  return(output)

}
