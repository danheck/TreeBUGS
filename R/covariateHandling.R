


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
    if(any(is.na(parSel)))
      stop("Check parameter names in covStructure. Problematic right now:\n  ", covTable$Parameter[is.na(parSel)])
    covTable$theta <- thetaNames$theta[parSel]
    covTable$Parameter <- thetaNames$Parameter[parSel]
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

######################## generate appropriate JAGS model string for latent-trait MPT
# covTable: result from preprocessing using covHandling()
# S: number of free parameters
covStringTrait <- function(covTable, S){
  modelString <- "\n## Probit Transformation and Covariate Handling ##\n"
  covPars <- c()


  ##### no predictors/ covariates: simple probit transformation as before
  if(is.null(covTable)){
    modelString <- "
for(i in 1:subjs) {
  for(s in 1:S){
    theta[s,i] <- phi(mu[s] + xi[s]*delta.part.raw[s,i])
  }
}
"

  ##### predictors included!
  }else{
    # which parameters include covariates?
    parWithCov <- unique(covTable$theta)
    modelString <- "\nfor(i in 1:subjs) {"
    for(s in 1:S){

      #### parameter WITH covariate
      if(s %in% parWithCov){

        selectLines <- covTable$theta == s
        # beginning of line as usual:
        modelString <- paste0(modelString,
                              "\ntheta[",s,",i] <- phi(mu[",s,"] + xi[",s,"]*delta.part.raw[",s,",i]")
        #
        for(cc in 1:sum(selectLines)){
          covIdx <-  covTable$covIdx[selectLines][cc]
          covParTmp <- paste0("slope_",
                                       covTable$Parameter[selectLines][1],"_",   # parameter label
                                       covTable$Covariate[selectLines][cc])      # covaraite label
          covPars <- c(covPars, covParTmp)
          modelString <- paste0(modelString, " + ",covParTmp,"*covData[i,",covIdx,"]")

        }
        modelString <- paste0(modelString, ")")

      }else{
        #### parameter WITHOUT covariate: as usual (but with explicit number as index)
        modelString <- paste0(modelString,
                              "\ntheta[",s,",i] <- phi(mu[",s,"] + xi[",s,"]*delta.part.raw[",s,",i])")
      }
    }
    modelString <- paste0(modelString, "\n}\n")

    # hyperpriors for slopes
    for(pp in 1:length(covPars)){
      modelString <- paste0(modelString, "\n", covPars[pp], " ~ dnorm(0,1)")
    }
  }


# lines in JAGS:
# additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
###################  ###################   ###################

return(list(modelString = modelString, covPars = covPars))
}
