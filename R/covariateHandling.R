


################### covariate handling of: covData, covStructure   ###################

# test files (Daniel):
# covData <- "../first_tests/covDatFac.csv"
# covData <- "../first_tests/covDat.csv"
# debug(TreeBUGS:::covHandling)



# covHandling: returns a nice table that assigns parameters to covariates (including indices)
#
# thetanames: data frame with assignment of parameters to indices / restrictions etc.
# covType: character vector with values "c"/"f"/"r" for continuous/fixed/random effects
covHandling <- function(covData, covStructure=NULL, N, thetaNames, covType=NULL){

  # get a clean list
  if(!is.null(covData)){

    # list / path to file
    if(is.character(covData)){
      covData <- read.csv(covData, header=T, sep= ";", strip.white = T)
    }
    if(is.character(covStructure)){
      covStructure <- as.list(read.csv(covStructure, header=F,stringsAsFactors=F, sep="}")$V1)
    }
    if(nrow(covData) != N){
      stop("Number of individuals in 'data' and 'covData' differs!")
    }

    covNames <- colnames(covData)
    if(is.null(covNames))
      stop("Check names of covariates in covData!")



    if(is.null(covStructure)){
      # default: all covariates for all parameters included
      covStructure <- lapply(unique(thetaNames$Parameter), function(tt) paste(tt, ";", paste(covNames, collapse=" ")))
    }

    # default: continuous / random covariates
    if(missing(covType) | is.null(covType)){
      cov.class <- sapply(covData, class)
      covType <- ifelse(cov.class %in% c("character", "factor"), "f", "c")
    }
    for(cc in 1:length(covType)){
      if(!all(covType %in% c("f","c","r"))  |  length(covType) != ncol(covData)){
        stop("Check definition of covType: should be a vector of the same length\n  ",
             "as there are columns in covData. Possible values are:\n  ",
             "'c' (continuous covariate), 'f' (fixed effect), or 'r' (random effect).")
      }
    }

    # recode factors from character values => integer with number of factor level as index
    covTmp <- covRecodeFactor(covData, covType)
    covData <- covTmp$covData
    covFactorLevels <- covTmp$covFactorLevels
    names(covFactorLevels) <- colnames(covData)

    # set up table and iterate across all combinations of parameters and covariates
    covTable <- data.frame()
    for(i in 1:length(covStructure)){
      sss <- strsplit(covStructure[[i]], ";")[[1]]
      pars <- strsplit(sss[1], " +")[[1]]
      covs <- strsplit(sss[2], " +")[[1]]
      pars <- pars[pars != ""]
      covs <- covs[covs != ""]

      # parameters:
      for(pp in 1:length(pars)){
        # covariates:
        for(cc in 1:length(covs)){
          covTable <- rbind(covTable, data.frame(Parameter = pars[pp], Covariate = covs[cc]))
        }
      }
    }

    # replace the constrained parameters by free parameters and remove redundant rows
    parSel <- match(covTable$Parameter,  thetaNames$Parameter)
    if(any(is.na(parSel)))
      stop("Check parameter names in covStructure. Problematic right now:\n  ",
           covTable$Parameter[is.na(parSel)])
    covTable$theta <- thetaNames$theta[parSel]
    covTable$Parameter <- thetaNames$Parameter[parSel]
    covTable$covIdx <- (1:ncol(covData))[match(covTable$Covariate,  covNames)]
    covTable <- covTable[!duplicated(paste(covTable$theta, covTable$covIdx)),]

    covTable$covType <- covType[covTable$covIdx]


    # get parameter labels for JAGS
    covTable$prefix <- ifelse(covTable$covType=="c", "slope", "factor")
    covTable$covPar <- apply(covTable[,c("prefix", "Parameter", "Covariate")],
                             1, paste, collapse="_")

  }else{
    # not covariates
    covTable <- NULL
    covFactorLevels <- NULL
  }

  return(list(covTable = covTable,
              covData = covData,
              covFactorLevels = covFactorLevels,
              covType = covType))
}




covStringBeta <- function(covTable){
  modelString <- "\n## Covariate Handling ##\n"
  # covTable$covPars <- NA


  if(!is.null(covTable)){
    pars <- unique(covTable$Parameter)
    modelString <- paste0(modelString,
           "for(i in 1:S){
  thetaSD[i] <- sd(theta[i,])
}
")
    cnt <- 0
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
        covTable$covPar[cnt<-cnt+1] <-  paste0("cor_", pars[pp],"_",covs[cc])
      }
    }
  }


  # lines in JAGS:
  # additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
  ###################  ###################   ###################

  return(list(modelString = modelString, covPars = covTable$covPar))
}



######################## generate appropriate JAGS model string for latent-trait MPT
# covTable: result from preprocessing using covHandling()
# S: number of free parameters
# covFactorLevels: list with factor levels for discrete covariates

covStringTrait <- function(covTable, S, covFactorLevels=NULL){
  modelString <- "\n## Probit Transformation and Covariate Handling ##\n"

  ##### no predictors/ covariates: simple probit transformation as before
  if(is.null(covTable)){
    covPars <- c() ; X_list <- list()
    modelString <- "
for(i in 1:subjs) {
  for(s in 1:S){
    theta[s,i] <- phi(mu[s] + xi[s]*delta.part.raw[s,i])
  }
}
"

  ############################### PREDICTORS IN PHI() TRANSFORMATION ########################
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

        # loop across covariates
        for(cc in 1:sum(selectLines)){
          covIdx <-  covTable$covIdx[selectLines][cc]

          if(covTable$covType[selectLines][cc] == "c"){
            ####### continuous covariate
            # add to model string:
            modelString <- paste0(modelString, " + ",covTable$covPar[selectLines][cc],
                                  "*covData[i,",covIdx,"]")

          }else{
            ####### discrete covariate
            modelString <- paste0(modelString, " + ",covTable$covPar[selectLines][cc],
                                  "[covData[i,",covIdx,"]]")
          }

        }
        modelString <- paste0(modelString, ")")

      }else{
        #### parameter WITHOUT covariate: as usual (but with explicit number in JAGS as index)
        modelString <- paste0(modelString,
                              "\ntheta[",s,",i] <- phi(mu[",s,"] + xi[",s,"]*delta.part.raw[",s,",i])")
      }
    }
    modelString <- paste0(modelString, "\n}\n")


    ############################### HYPERPRIORS     ###############################

    covPars <- covTable$covPar
    X_list <- list()  # list with design matrices for fixed effects)
    for(pp in 1:nrow(covTable)){
      if(covTable$covType[pp] == "c"){

        # hyperpriors for slopes
        modelString <- paste0(modelString, "\n", covPars[pp], " ~ dnorm(0,1)")

      }else if(covTable$covType[pp] == "r"){
        numLevel <- length(covFactorLevels[[covTable$covIdx[pp]]])
        # random effects:  (cf. Rouder et al (2012))
        modelString <- paste0(modelString,
                              "\nfor(level in 1:", numLevel,"){\n  ",
                              covPars[pp], "[level] ~ dnorm(0, tau_",covPars[pp],") \n}\n",
                              "tau_", covPars[pp], " ~ dchisq(1)\n",
                              "SD_", covPars[pp], " <- sqrt(inverse(tau_",covPars[pp] ,"))")

        covPars <- c(covPars, paste0("SD_", covPars[pp]))
      }else{
        # fixed effects: (Rouder et al, 2012)
        numLevel <- length(covFactorLevels[[covTable$covIdx[pp]]])
        Z <- diag(numLevel) - 1/numLevel

        # add design matrix to list of design matrices:
        X_list <- c(X_list,
                    list(X = eigen(Z, symmetric=T)$vectors[,1:(numLevel-1)]))
        names(X_list)[length(X_list)] <- paste0("X_", covPars[pp])

        # model string: only (numLevel-1) priors
        modelString <- paste0(modelString,

                              # univariate normal priros on numLevel-1 of stransformed parameters
                              "\nfor(level in 1:", numLevel-1, "){\n  ",
                              "s",covPars[pp],"[level] ~ dnorm(0, tau_",covPars[pp],") \n}",

                              # contrast coding of fixed effects:
                              "\nfor(level in 1:", numLevel, "){\n  ",
                              covPars[pp],"[level] <- inprod(",
                              "s",covPars[pp],", X_",covPars[pp],"[level,])\n} \n",

                              # hyperpriors:
                              "tau_", covPars[pp], " ~ dchisq(1)\n",
                              "SD_", covPars[pp], " <- sqrt(inverse(tau_", covPars[pp],"))")
        # before calling sampler:        assign(paste0("X_",names(X_list)[pp]), X_list[[pp]])

        covPars <- c(covPars, paste0("SD_", covPars[pp]))
      }
    }
  }


# lines in JAGS:
# additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
###################  ###################   ###################

return(list(modelString = modelString, covPars =  covPars, X_list=X_list))
}
