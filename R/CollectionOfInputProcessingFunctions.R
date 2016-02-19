# #' Unique branches for each Tree with the same answer and sum the corresponing formulas
# #'
# #' @param Tree Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
mergeBranches<-function(Tree){ # OLD ,DataNames){ # Unique branches for each Tree with one answer and sum the corresponing formulas

  treeNames <- sort(unique(Tree$Tree))
  catNames <- as.list(by(Tree$Category, Tree$Tree, function(xx) sort(unique(xx))))
  names(catNames) <- treeNames

  NewTree <- data.frame(Tree = rep(treeNames, sapply(catNames, length)),
                        Category = as.character(unlist(catNames)),
                        Equation = "...", stringsAsFactors = FALSE)
  for(tt in 1:length(treeNames)){
    for(cc in 1:length(catNames[[tt]])){
      selTree <- Tree$Tree == treeNames[tt] & Tree$Category == catNames[[tt]][cc]
      selNewTree <- NewTree$Tree == treeNames[tt] & NewTree$Category == catNames[[tt]][cc]
      NewTree$Equation[selNewTree] <- paste(Tree[selTree,]$Equation, collapse="+")
    }
  }

  return(NewTree)
}



# #' Extract all parameter from the formulas of a given model
# #'
# #' @param TreeData Data returned from readEQN()
# #' @author Nina R. Arnold, Denis Arnold
# #' @export
getParameter<-function(TreeData){

  Parameter <- unique(unlist(strsplit(TreeData$Equation,
                                      split="\\*|\\(|\\)|\\-|\\+")))
  r <- c(which(nchar(Parameter)==0), grep("^[0-9]+$|^[0-9]+\\.[0-9]+", Parameter))
  Parameter <- Parameter[-r]
  #   Parameter=c(sort(Parameter[grepl("[A-Z]",Parameter)]),
  #               sort(Parameter[!grepl("[A-Z]",Parameter)]))
  # returns errors if model does not contain uppercase parameters

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
      stop("Number of categories differs in eqn and csv file.")
    }else{
      stop("At least one name of the categories differs in eqn and csv file")
    }
  }

  data <- data[,Category] #order data columns according to Tree category label order

  return(data)
}

