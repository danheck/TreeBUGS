


################### covariate handling of: covData, covStructure   ###################
# covData <- "../first_tests/covDat.csv"
# debug(TreeBUGS:::covHandling)
covHandling <- function(covData, covStructure=NULL, N, thetaNames){

  # get a clean list
  if(!is.null(covData)){

    # list / path to file
    if(is.character(covData)){
      covData <- read.csv(covData, header=T, sep= ",",  strip.white = TRUE)
    }
    if(is.character(covStructure)){
      covStructure <- as.list(read.csv(covStructure, header=F,stringsAsFactors=F, sep="}")$V1)
    }
    if(nrow(covData) != N){
      stop("Number of individuals in 'data' and 'covData' differs!")
    }

    covNames <- colnames(covData)
    if(is.null(covNames))
      stop("Check names of covariates")



    if(is.null(covStructure)){
      # default: all covariates for all parameters included
      covStructure <- lapply(unique(thetaNames$Parameter), function(tt) paste(tt, ";", paste(covNames, collapse=" ")))
    }
    covTable <- data.frame()
    for(i in 1:length(covStructure)){
      sss <- strsplit(covStructure[[i]], ";")[[1]]
      pars <- strsplit(sss[1], " +")[[1]]
      covs <- strsplit(sss[2], " +")[[1]]
      pars <- pars[pars != ""]
      covs <- covs[covs != ""]
      for(pp in 1:length(pars)){
        for(cc in 1:length(covs)){
          covTable <- rbind(covTable, data.frame(Parameter = pars[pp], Covariate = covs[cc]))
        }
      }
    }

    # replace by constrained parameters and remove redundant rows
    parSel <- match(covTable$Parameter,  thetaNames$Parameter)
    covTable$theta <- thetaNames$theta[parSel]
    covTable$Parameter <- thetaNames$Parameter[covTable$theta]
    covTable$covIdx <- (1:ncol(covData))[match(covTable$Covariate,  covNames)]
    covTable <- covTable[!duplicated(paste(covTable$theta, covTable$covIdx)),]
  }else{
    covTable <- NULL
  }


  return(covTable)
}




covStringBeta <- function(covTable){
  modelString <- "\n## Covariate Handling ##\n"
  covPars <- c()


  if(!is.null(covTable)){
    pars <- unique(covTable$Parameter)
    modelString <- paste0(modelString,
           "for(i in 1:S){
  thetaSD[i] <- sd(theta[i,])
}
")
    for(pp in 1:length(pars)){
      sel <- covTable$Parameter == pars[pp]
      thetaIdx <- covTable$theta[sel][1]
      covs <- covTable$Covariate[sel]
      for(cc in 1:length(covs)){
        covIdx <-  covTable$covIdx[sel][cc]

        ### requires commputation of correlation using only sum / sd
        #  (par-mean(par))*(cov-mean(cov))/(sd(par) * sd(cov))

        modelString <- paste0(modelString,
                              "cor_", pars[pp],"_",covs[cc],
                              " <- mean( (theta[",thetaIdx, ",]-mean(theta[",thetaIdx,
                              ",]))*(covData[,",covIdx, "]-mean(covData[,",covIdx,
                              "])) )/covSD[",covIdx,
                              "] / thetaSD[",thetaIdx,"]\n")
# nice, not supported: cor(theta[",thetaIdx, ",], covData[,",covIdx, "])
        covPars <- c(covPars, paste0("cor_", pars[pp],"_",covs[cc]))
      }
    }
  }


  # lines in JAGS:
  # additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
  ###################  ###################   ###################

  return(list(modelString = modelString, covPars = covPars))
}
