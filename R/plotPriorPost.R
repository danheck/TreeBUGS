
#' Plot Prior vs. Posterior Distribution
#'
#' Allows to judge how much the data informed the parameter posterior distributions compared to the prior.
#' @inheritParams plotFit
#' @inheritParams priorPredictive
#' @inheritParams plotPrior
#' @param M number of random samples to approximate prior distributions
# @param scale whether to scale the plots to give a complete picture of the \code{"prior"} or posterior (\code{scale="post"})
#' @param ci credibility interval indicated by vertical red lines
#' @details Prior distributions are shown as blue, dashed lines, whereas posterior distributions are shown as solid, black lines.
#' @export
plotPriorPost <- function(fittedModel, probitInverse = "mean", M=2e5, ci=.95, nCPU=3,...){

  mfrow <- par()$mfrow

  S <- length(fittedModel$mptInfo$thetaUnique)
  samples <- sampleHyperprior(fittedModel$mptInfo$hyperprior, M=M, S=S,
                              probitInverse=probitInverse, nCPU=nCPU)

  for(s in 1:S){
    label <- ifelse( S== 1, "", paste0("[",s,"]"))
    mean.post <- unlist(fittedModel$runjags$mcmc[,paste0("mean", label)])
    d.mean <- density(mean.post, from=0, to=1, na.rm = TRUE)
    prior.mean <- density(samples$mean[,s], from=0,to=1, na.rm = TRUE)
    ci.mean <- quantile(mean.post, c((1-ci)/2,1-(1-ci)/2))

    if(fittedModel$mptInfo$model == "betaMPT"){
      sd.post <- unlist(fittedModel$runjags$mcmc[,paste0("sd", label)])
      ci.sd <- quantile(sd.post, c((1-ci)/2,1-(1-ci)/2))
      d.sd <- density(sd.post, from=0, to=.5, na.rm = TRUE)
      prior.sd <- density(samples$sd[,s], from=0, to=.5, na.rm = TRUE)
      xlab.sd = "Group SD (probability)"

      ####### traitMPT
    }else{
      sig.post <- unlist(fittedModel$runjags$mcmc[,paste0("sigma", label)])
      if(probitInverse == "mean_sd"){
        mean_sd <- probitInverse(qnorm(mean.post), sig.post)
        d.mean <- density(mean_sd[,"mean"], from=0, to=1, na.rm = TRUE)
        d.sd <- density(mean_sd[,"sd"], from = 0, to = .5, na.rm = TRUE)
        ci.sd <- quantile(mean_sd[,"sd"], c((1-ci)/2,1-(1-ci)/2))
        ci.mean <- quantile(mean_sd[,"mean"], c((1-ci)/2,1-(1-ci)/2))
        prior.sd <- density(samples$sd[,s], from=0, to=.5 ,na.rm = TRUE)
      }else{
        prior.sd <- density(samples$sd[,s], from = 0, na.rm = TRUE)
        d.sd <- density(sig.post, from = 0, na.rm = TRUE)
        ci.sd <- quantile(sig.post, c((1-ci)/2,1-(1-ci)/2))
        if(probitInverse == "none"){
          d.mean <- density(qnorm(mean.post), na.rm = TRUE)
          ci.mean <- quantile(qnorm(mean.post), c((1-ci)/2,1-(1-ci)/2))
          prior.mean <- density(samples$mean[,s], na.rm = TRUE)
        }
      }

      xlab.sd = ifelse(probitInverse=="mean_sd",
                       "Group SD (probability scale)",
                       "Group SD (probit scale)")
    }

    par(mfrow=1:2)
    tmp <- readline(prompt = "Press <Enter> to show the next plot.")
    plot(d.mean, main=paste0( "Group mean of ", fittedModel$mptInfo$thetaUnique[s]),
         xlab="Group mean", las=1, ...)
    lines(prior.mean, col="blue", lty="dashed")
    abline(v= ci.mean, col="red")
    plot(d.sd,   main=paste0("Group SD of ", fittedModel$mptInfo$thetaUnique[s]),
         xlab=xlab.sd, las=1, ...)
    lines(prior.sd, col="blue", lty="dashed")
    abline(v=ci.sd,col="red")
  }

  if(fittedModel$mptInfo$model == "traitMPT" & S>1){
    cnt <- 0
    for(s1 in 1:(S-1)){
      for(s2 in (s1+1):S){
        d.cor <- density(unlist(fittedModel$runjags$mcmc[,paste0("rho[",s1,",",s2,"]")]),
                         from=-1, to=1)
        prior.cor <- density(samples$rho[s1,s2,], from=-1, to=1)
        if(cnt/2 == round(cnt/2))
          tmp <- readline(prompt = "Press <Enter> to show the next plot.")
        cnt <- cnt+1
        plot(d.cor, main=paste0( "Correlation between ",
                                  fittedModel$mptInfo$thetaUnique[s1], " and ",
                                  fittedModel$mptInfo$thetaUnique[s2]),
             xlab="Correlation (on latent probit scale)", las=1, ...)
        lines(prior.cor, col="blue", lty="dashed")
        abline(v= quantile(unlist(fittedModel$runjags$mcmc[,paste0("rho[",s1,",",s2,"]")]),
                           c((1-ci)/2,1-(1-ci)/2)), col="red")
      }
    }
  }

  par(mfrow=mfrow)
}

