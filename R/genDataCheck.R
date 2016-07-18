
################### HELPER FUNCTIONS ####################################


checkNumItems <- function(numItems, treeLabels){
  if(any(sort(treeLabels) != sort(names(numItems))))
    stop("Names for numItems do not match the tree labels in EQN file:\n  ",
         paste(substr(treeLabels, 3, 100), collapse=", "))
  if(length(numItems) != length(treeLabels))
    stop("Argument numItems has the wrong length (should be", length(treeLabels),")")
}

checkThetaNames <- function(theta, thetaNames){
  if(any(sort(thetaNames) != sort(colnames(theta))))
    stop("Column names of theta do not match parameters in EQN file:\n  ",
         paste(thetaNames, collapse=", "))
}

checkNamingMatrix <- function(S, thetaNames, matrix, matrixName = "rho"){
  if(any(S != dim(matrix)))
    stop("Dimensions of matrix '",matrixName,"' not correct, should be ", S, "x", S)

  if(is.null(dimnames(matrix))){
    warning("Matrix '",matrixName,"' not named. Internal order of parameters is used.\n  See ?readMultiTree and check parameters by generatedData$parameters")
    dimnames(matrix) <- list(thetaNames, thetaNames)
  }else if(any(sort(thetaNames) != sort(rownames(matrix)))){
    stop("Row names of matrix '", matrixName,"' do not match parameter labels in eqn file.")
  }else if(any(sort(thetaNames) != sort(colnames(matrix)))){
    stop("Column names of matrix '", matrixName,"' do not match parameter labels in eqn file.")
  }else{
    matrix <- matrix[thetaNames,]
    matrix <- matrix[,thetaNames]
  }

  if(any(diag(matrix) != 1))
    stop("Diagonal must contain ones!")

  if(any(matrix != t(matrix)))
    stop("Matrix ", matrixName, " must be symmetric!")

  if(any(matrix < -1 | matrix >1))
    stop("'",matrixName,"' cannot be negative!")

  return(matrix)
}

checkNaming <- function(S, thetaNames, vector, vectorName){
  if(S != length(vector))
    stop("Length of '",vectorName,"' not correct, should be ", S)
  if(is.null(names(vector))){
    warning("Vector '",vectorName,"' not named. Internal order of parameters is used, see ?readMultiTree and check parameters by generatedData$parameters")
    names(vector) <- thetaNames
    #     cat(vectorName, ":\n")
    #     print(vector)
  }else if(any(sort(thetaNames) != sort(names(vector)))){
    stop("Parameter names of vector '",vectorName,"' do not match parameter labels in eqn file.")
  }else{
    vector <- vector[thetaNames]
  }

  if(any(vector<0))
    stop("'",vectorName,"' cannot be negative!")

  return(vector)
}
