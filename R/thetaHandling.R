
# #'  Set constants, replace parameters for thetas
# #'
# #' @param mergedTree is the output of mergeBranches
# #' @param restrictions Either the (full path to the) file that specifies which parameters should be constants and which should # be equal; or a list of parameter restrictions
# #' @author Nina R. Arnold, Denis Arnold, Daniel Heck
# #' @export
thetaHandling <-function(mergedTree, restrictions){

  mergedTree$Eqn <- mergedTree$Equation
  Parameter <- getParameter(mergedTree)

  SubPar <- data.frame(Parameter = Parameter,
                       theta = 1:length(Parameter),
                       sub="", stringsAsFactors=FALSE)

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
    restrVector <- gsub(" ", "", restrVector, fixed = TRUE)

    for(k in 1:length(restrVector)){
      splitRestr <- strsplit(restrVector[k], "=")[[1]]
      selFE <- splitRestr %in% c("FE", "FIXED_EFFECT")
      if(length(splitRestr) == 1){
        warning("Restriction not well defined: Equality sign '=' missing in:\n  ",splitRestr)

      }else if(any(selFE)){
        # fixed effect: negative index
        index <- - match(splitRestr[!selFE], SubPar$Parameter)
        SubPar[ - index, "theta"]<- index[1] - 1 # substract one to make it different from equality constraints

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

  isConstant <- SubPar$theta <= 0 & SubPar$theta >= -1
  isFE <-  SubPar$theta < -1

  ############## replaced constant parameter values:
  if(any(isConstant)){
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
  }else{
    constants <- NULL
  }
  if(any(isFE)){
    fixedPar <- SubPar[isFE,, drop=FALSE]
    # fixedPar$theta <- - fixedPar$theta - 1
    SubPar <- SubPar[!isFE,, drop=FALSE]
    tmp <- unique(fixedPar$theta)
    for(tt in 1:length(tmp)){
      fixedPar$theta[fixedPar$theta == tmp[tt]] <- tt
    }
    fixedPar$sub <- paste0("thetaFE[", fixedPar$theta , "]")


    for(j in 1:nrow(fixedPar)){
      mergedTree$Equation <- gsub(
        pattern = paste0("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
                         fixedPar$Parameter[j],
                         "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
        replacement = paste0("\\1",fixedPar$sub[j],"\\2"),
        mergedTree$Equation, perl=T)
    }
  } else{
    fixedPar <- NULL
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

    # two times for special case that same parameter occurs twice directly after each other: u*u
    mergedTree$Equation <- gsub(
      pattern = paste0("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",
                       SubPar$Parameter[j],
                       "(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)"),
      replacement = paste0("\\1",SubPar$sub[j],"\\2"),
      mergedTree$Equation, perl=T)
  }

  output=list(SubPar = SubPar, mergedTree = mergedTree,
              constants=constants, fixedPar=fixedPar,
              restrictions = restrictions)

  return(output)

}
