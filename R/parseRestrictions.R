

parseRestrictions <- function(mpt, restrictions){

  # thetaNames <- data.frame(Parameter = colnames(mpt$a),
  #                          theta=1:nrow(mpt$a))

  mpt$c <- rep(1, nrow(mpt$a))

  ############################### constraints in MPT file (e.g., .5*Do  )

  suppressWarnings(parConst <- as.numeric(colnames(mpt$a)))
  numConst <- sum(!is.na(parConst))
  if(numConst>0){
    idx <- which(!is.na(parConst))
    for(s in 1:numConst){
      mpt$c <- mpt$c *
        parConst[idx[s]]^mpt$a[,idx[s]] *
        (1-parConst[idx[s]])^mpt$b[,idx[s]]
    }
    mpt$a <- mpt$a[,-idx,drop=FALSE]
    mpt$b <- mpt$b[,-idx,drop=FALSE]
  }

  ############################### Contraints in list "restrictions"

  parLabels <- colnames(mpt$a)

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
      if(length(splitRestr) == 1){
        warning("Restriction not well defined: Equality sign '=' missing in:\n  ",splitRestr)

      }else{
        ######### equality constraints

        index <- match(splitRestr, colnames(mpt$a))
        suppressWarnings(consts <- as.numeric(splitRestr))

        if(  all(is.na(consts)) ){
          # only parameters without constants
          if(any(is.na(index))){
            error <- paste0("Restriction contains parameters not contained in the model:\n  ",
                            paste(splitRestr, collapse="="))
            stop(error)
          }

          # thetaNames$theta[index[2:length(index)]] <- index[1]

          # replace index
          # mpt$a[,index[1]] <- rowSums(mpt$a[,index])
          # mpt$b[,index[1]] <- rowSums(mpt$b[,index])
          # mpt$a <- mpt$a[,-index[2:length(index)],drop =FALSE]
          # mpt$b <- mpt$b[,-index[2:length(index)],drop =FALSE]
          mpt$a[,min(index)] <- rowSums(mpt$a[,index])
          mpt$b[,min(index)] <- rowSums(mpt$b[,index])
          mpt$a <- mpt$a[,-index[-which.min(index)],drop =FALSE]
          mpt$b <- mpt$b[,-index[-which.min(index)],drop =FALSE]


        }else if(sum(!is.na(consts)) == 1){
          # contrained to a single constant value
          CONST <- consts[!is.na(consts)]
          if(CONST <0 | CONST >1){
            error <- paste0("Check parameter restrictions. Constants are not in the interval [0,1]: ",
                            restrVector[k])
            warning(error)
          }

          mpt$c <- mpt$c *
            apply(CONST^mpt$a[,index[!is.na(index)],drop=FALSE] *
                    (1-CONST)^mpt$b[,index[!is.na(index)],drop=FALSE], 1, prod)

          mpt$a <- mpt$a[,-index[!is.na(index)],drop =FALSE]
          mpt$b <- mpt$b[,-index[!is.na(index)],drop =FALSE]

          # thetaNames <- thetaNames[,]  $theta[index[2:length(index)]] <- index[1]

        }else{
          stop("Restrictions should not contain more than one constant!")
        }
      }
    }
  }
  # mpt$thetaNames <- thetaNames

  mpt
}
