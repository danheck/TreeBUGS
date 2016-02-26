
# make index variables to compute T1 statistic by group
getGroupT1 <- function(covData, predType, T1group=NULL){

  if(!is.null(T1group)){
    if(!T1group %in% colnames(covData)){
      stop("T1group not found in covData (matching by column names)")
    }
    if(is.null(covData) || is.null(T1group)){
      return(NULL)
    }

    groupTab <- table(covData[,T1group])
    NgroupT1 <- as.vector(groupTab)
    groupNames <- names(groupTab)
    names(NgroupT1) <- groupNames
    N <- sum(NgroupT1)
    G <- length(NgroupT1)

    #### general index vector, stored in matrix
    groupMatT1 <- matrix(0, G, max(NgroupT1),
                         dimnames=list(groupNames, NULL))
    cnt <- rep(0, G)
    for(i in 1:N){
      group <- covData[i,T1group]
      groupMatT1[group,cnt[group] <- cnt[group]+1] <- i
    }
    list(groupMatT1=groupMatT1,
         NgroupT1=NgroupT1)
  }else{
    NULL
  }
}
