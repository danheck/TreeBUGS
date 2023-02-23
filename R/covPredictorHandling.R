

######################## generate appropriate JAGS model string for latent-trait MPT
# covTable: result from preprocessing using covHandling()
# S: number of free parameters
# predFactorLevels: list with factor levels for discrete covariates

covStringPredictor <- function(covTable, S, predFactorLevels = NULL, IVprec = "dgamma(.5,.5)") {
  modelString <- "\n## Probit Transformation and Covariate Handling ##\n"

  ##### no predictors/ covariates: simple probit transformation as before
  if (is.null(covTable)) {
    covPars <- c()
    X_list <- list()
    modelString <- "
    for(i in 1:subjs) {
      for(s in 1:S){
        theta[s,i] <- phi(mu[s] + xi[s]*delta.part.raw[s,i])
      }
    }
    "

    ############################### PREDICTORS IN PHI() TRANSFORMATION ########################
  } else {
    # which parameters include covariates?
    parWithCov <- unique(covTable$theta)
    modelString <- "\nfor(i in 1:subjs) {"
    for (s in 1:S) {
      #### parameter WITH covariate
      if (s %in% parWithCov) {
        selectLines <- covTable$theta == s
        # beginning of line as usual:
        modelString <- paste0(
          modelString,
          "\ntheta[", s, ",i] <- phi(mu[", s, "] + xi[", s, "]*delta.part.raw[", s, ",i]"
        )

        # loop across covariates
        for (cc in 1:sum(selectLines)) {
          covIdx <- covTable$covIdx[selectLines][cc]

          if (covTable$predType[selectLines][cc] == "c") {
            ####### continuous covariate
            # add to model string:
            modelString <- paste0(
              modelString, " + ", covTable$covPar[selectLines][cc],
              "*covData[i,", covIdx, "]"
            )
          } else {
            ####### discrete covariate
            modelString <- paste0(
              modelString, " + ", covTable$covPar[selectLines][cc],
              "[covData[i,", covIdx, "]]"
            )
          }
        }
        modelString <- paste0(modelString, ")")
      } else {
        #### parameter WITHOUT covariate: as usual (but with explicit number in JAGS as index)
        modelString <- paste0(
          modelString,
          "\ntheta[", s, ",i] <- phi(mu[", s, "] + xi[", s, "]*delta.part.raw[", s, ",i])"
        )
      }
    }
    modelString <- paste0(modelString, "\n}\n")

    # parameter labels in JAGS:
    covPars <- covTable$covPar



    ############################### standardization of regression slopes (--> to get correlations)

    if (!is.null(covTable)) {
      modelString <- paste0(
        modelString,
        "\n### standardization of regression slopes",
        "(only for continuous predictors)\n\n"
      )
      select_cont <- which(covTable$predType == "c")
      slope_std <- gsub("slope_", "slope_std_", covTable$covPar[select_cont], fixed = TRUE)
      covPars <- c(covPars, slope_std)
      for (pp in seq_along(slope_std)) {
        s <- covTable$theta[select_cont[pp]]
        k <- covTable$covIdx[select_cont[pp]]

        # for partially standardized slopes (z-standardized predictors):
        # slope_std[pp] <- slope[pp] / sqrt(slope[pp]^2 + sigma[s]^2)

        # for unstandardized predictors:
        # slope_std[pp] <- slope[pp] * sqrt(covVar[k]) / sqrt(slope[pp]^2 + covVar[k] + sigma[s]^2)

        modelString <- paste0(
          modelString,
          slope_std[pp], " <- ",
          covTable$covPar[select_cont[pp]], " * sqrt(covVar[", k, "]) / ",
          "sqrt(",
          "pow(", covTable$covPar[select_cont[pp]], ",2)",
          "* covVar[", k, "]",
          " + Sigma[", s, ",", s, "]",
          ") \n"
        )
      }
      modelString <- paste0(modelString, "\n")
    }


    ############################### HYPERPRIORS     ###############################

    X_list <- list() # list with design matrices for fixed effects)

    # inverse-g prior for Zellner-Siow priors (contin. predictors)
    if ("c" %in% covTable$predType) {
      parWithPred <- unique(subset(covTable, covTable$predType == "c")$Parameter)
      for (pp in seq_along(parWithPred)) {
        modelString <- paste0(
          modelString, "\n",
          "ginv_", parWithPred[pp],
          ifelse(is.numeric(IVprec), " <- ", " ~ "), IVprec
        )
      }
    }

    for (pp in 1:nrow(covTable)) {
      if (covTable$predType[pp] == "c") {
        # hyperpriors for slopes
        modelString <- paste0(
          modelString, "\n",
          # "tau_",covPars[pp]," <- ginv_",covTable$Parameter[pp],
          # "* covVar[",covTable$covIdx[pp],"]\n",
          # * sigma[",covTable$theta[pp],"]^(-2)

          covPars[pp], " ~ dnorm(0,", "ginv_",
          covTable$Parameter[pp], "* covVar[", covTable$covIdx[pp], "])"
        )
      } else if (covTable$predType[pp] == "r") {
        numLevel <- length(predFactorLevels[[covTable$covIdx[pp]]])
        # random effects:  (cf. Rouder et al (2012))
        modelString <- paste0(
          modelString,
          "\nfor(level in 1:", numLevel, "){\n  ",
          covPars[pp], "[level] ~ dnorm(0, tau_", covPars[pp], ") \n}\n",
          "tau_", covPars[pp], " ~ dgamma(.5, .5)\n",
          "SD_", covPars[pp], " <- sqrt(inverse(tau_", covPars[pp], "))"
        )

        covPars <- c(covPars, paste0("SD_", covPars[pp]))
      } else {
        # fixed effects: (Rouder et al, 2012)
        numLevel <- length(predFactorLevels[[covTable$covIdx[pp]]])
        Z <- diag(numLevel) - 1 / numLevel

        # add design matrix to list of design matrices:
        X_list <- c(
          X_list,
          list(X = as.matrix(eigen(Z, symmetric = T)$vectors[, 1:(numLevel - 1)]))
        )
        names(X_list)[length(X_list)] <- paste0("X_", covPars[pp])

        # model string: only (numLevel-1) priors
        modelString <- paste0(
          modelString,

          # univariate normal priros on numLevel-1 of stransformed parameters
          "\nfor(level in 1:", numLevel - 1, "){\n  ",
          "s", covPars[pp], "[level] ~ dnorm(0, tau_", covPars[pp], ") \n}",

          # contrast coding of fixed effects:
          "\nfor(level in 1:", numLevel, "){\n  ",
          covPars[pp], "[level] <- inprod(",
          "s", covPars[pp], ", X_", covPars[pp], "[level,])\n} \n",

          # hyperpriors:
          "tau_", covPars[pp], " ~ dgamma(.5, .5)\n",
          "SD_", covPars[pp], " <- sqrt(inverse(tau_", covPars[pp], "))"
        )
        # before calling sampler:        assign(paste0("X_",names(X_list)[pp]), X_list[[pp]])

        covPars <- c(covPars, paste0("SD_", covPars[pp]))
      }
    }
  }


  # lines in JAGS:
  # additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
  ###################  ###################   ###################

  list(modelString = modelString, covPars = covPars, X_list = X_list)
}
