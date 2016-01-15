
# factor handling

# recodes factor labels into integer index
# returns a pure numeric data set covData and the old factor labels
covRecodeFactor <- function(covData, covType){

  numCov <- ncol(covData)
  covFactorLevels <- vector("list", numCov)

  for(i in 1:numCov){

    if(covType[i] != "c"){

      # store factor level labels:
      covFactorLevels[[i]] <- sort(unique(levels(as.factor(covData[,i]))))
      nLevel <- length(covFactorLevels[[i]])

      if(nLevel <= 1){
        stop("Factor", colnames(covData)[i], "has only one factor level!")
      }

      # replace factor levels by an integer index:
      covData[,i] <- match( as.factor(covData[,i]), covFactorLevels[[i]] )
    }
  }


  list(covData = covData, covFactorLevels = covFactorLevels)
}

