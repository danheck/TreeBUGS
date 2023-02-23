################### HELPER FUNCTIONS ####################################

checkNumItems <- function(numItems, treeLabels) {
  if (any(sort(treeLabels) != sort(names(numItems)))) {
    stop(
      "Names for numItems do not match the tree labels in EQN file:\n  ",
      paste(substr(treeLabels, 3, 100), collapse = ", ")
    )
  }
  if (length(numItems) != length(treeLabels)) {
    stop("Argument numItems has the wrong length (should be", length(treeLabels), ")")
  }
  numItems[treeLabels]
}

checkThetaNames <- function(theta, thetaNames) {
  if (any(sort(thetaNames) != sort(colnames(theta)))) {
    stop(
      "Column names of theta do not match parameters in EQN file:\n  ",
      paste(thetaNames, collapse = ", ")
    )
  }
  theta[, thetaNames, drop = FALSE]
}

checkNamingMatrix <- function(S, thetaNames, matrix,
                              matrixName = "rho", warning = TRUE) {
  if (any(S != dim(matrix))) {
    stop("Dimensions of matrix '", matrixName, "' not correct, should be ", S, "x", S)
  }

  if (is.null(dimnames(matrix))) {
    if (warning) {
      warning(
        "Matrix '", matrixName, "' not named. Internal order of parameters is used.\n",
        "See ?readMultiTree and check parameters by generatedData$parameters"
      )
    }
    dimnames(matrix) <- list(thetaNames, thetaNames)
  } else if (any(sort(thetaNames) != sort(rownames(matrix)))) {
    stop("Row names of matrix '", matrixName, "' do not match parameter labels in eqn file.")
  } else if (any(sort(thetaNames) != sort(colnames(matrix)))) {
    stop("Column names of matrix '", matrixName, "' do not match parameter labels in eqn file.")
  } else {
    matrix <- matrix[thetaNames, , drop = FALSE]
    matrix <- matrix[, thetaNames, drop = FALSE]
  }

  if (any(diag(matrix) != 1)) {
    stop("Diagonal must have ones!")
  }

  if (any(abs(matrix - t(matrix)) > 1e-10)) {
    stop("Matrix ", matrixName, " must be symmetric!")
  }

  if (any(matrix < -1 | matrix > 1)) {
    stop("'", matrixName, "' cannot be negative!")
  }

  return(matrix)
}

checkNaming <- function(S, thetaNames, vector, vectorName,
                        interval = c(0, Inf), warning = TRUE) {
  if (S != length(vector)) {
    stop("Length of '", vectorName, "' not correct, should be ", S)
  }
  if (is.null(names(vector))) {
    if (warning) {
      warning(
        "Vector '", vectorName, "' not named. Internal order of parameters",
        " is used, see ?readMultiTree and check parameters by generatedData$parameters"
      )
    }
    names(vector) <- thetaNames
  } else if (any(sort(thetaNames) != sort(names(vector)))) {
    stop("Parameter names of vector '", vectorName, "' do not match parameter labels in eqn file.")
  } else {
    vector <- vector[thetaNames]
  }

  if (any(vector < interval[1] | vector > interval[2])) {
    stop("'", vectorName, "' cannot be below ", interval[1], " or above ", interval[2], ".")
  }

  return(vector)
}
