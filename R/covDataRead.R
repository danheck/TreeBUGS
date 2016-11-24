
# read, check and mean-center covData
covDataRead <- function(covData, N, binaryToNumeric=FALSE){

  if(!is.null(covData)){

    # list / path to file
    if(is.character(covData)){
      covData <- read.csv(covData, header=T, sep= ",", strip.white = T)
    }
    try(covData <- as.data.frame(covData))


    if(nrow(covData) != N){
      stop("Number of individuals in 'data' and 'covData' differs!")
    }

    if(is.null(colnames(covData)))
      stop("Check names of covariates in covData!")

    if(binaryToNumeric){
      for(k in 1:ncol(covData)){
        if(class(covData[,k]) %in% c("factor","ordered","character")){
          if(length(unique(covData[,k])) == 2){
            covData[,k] <- as.numeric(as.factor(covData[,k]))
          }
        }
      }
    }
  }
  covData
}



# get default values for predType
predTypeDefault <- function(covData, predType=NULL){
  if(!is.null(covData)){
    # default: continuous / random covariates
    if(missing(predType) || is.null(predType)){
      cov.class <- sapply(covData, class)
      predType <- ifelse(cov.class %in% c("character", "factor"), "f",
                        ifelse(cov.class %in% c("integer", "numeric","matrix"), "c",""))
    }
    for(cc in 1:length(predType)){
      if(!all(predType %in% c("f","c","r"))){
        stop("Check definition of predType: should be a vector of the same length\n  ",
             "as there are columns in covData. Possible values are:\n",
             "  'c' (continuous variable),",
             "\n  'f' (fixed effects factor; only traitMPT), or",
             "\n  'r' (random effects factor; only traitMPT).")
      }
    }
  }else{
    predType <- NULL
  }
  predType
}


# mean-centered variables as default (does not matter much for correlational analyses)
covDataCenter <- function(covData, predType){
  if(!is.null(covData)){
    for(i in 1:ncol(covData)){
      if(predType[i] == "c"){
        scaled <- scale(covData[,i], center=TRUE, scale = FALSE)  # centering of continuous variables
        if(any(scaled != covData[,i])){
          covData[,i] <- c(scaled)
        }
      }
    }
  }

  covData
}
