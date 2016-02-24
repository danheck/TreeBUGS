
# factor handling

# recodes factor labels into integer index
# returns a pure numeric data set covData and the old factor labels
covRecodeFactor <- function(covData, predType){

  numCov <- ncol(covData)
  predFactorLevels <- vector("list", numCov)

  for(i in 1:numCov){

    if(predType[i] != "c"){

      # store factor level labels:
      predFactorLevels[[i]] <- sort(unique(levels(as.factor(covData[,i]))))
      nLevel <- length(predFactorLevels[[i]])

      if(nLevel <= 1){
        stop("Factor", colnames(covData)[i], "has only one factor level!")
      }

      # replace factor levels by an integer index:
      covData[,i] <- match( as.factor(covData[,i]), predFactorLevels[[i]] )
    }
  }


  list(covData = covData, predFactorLevels = predFactorLevels)
}

