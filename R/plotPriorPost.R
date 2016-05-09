
#' Plot Prior vs. Posterior Distribution
#'
#' Allows to judge how much the data informed the parameter posterior distributions compared to the prior.
#' @inheritParams plotFit
#' @param M number of random samples to approximate prior distributions
# @param scale whether to scale the plots to give a complete picture of the \code{"prior"} or posterior (\code{scale="post"})
#' @param ci credibility interval indicated by vertical red lines
#' @details Prior distributions are shown as blue, dashed lines, whereas posterior distributions are shown as solid, black lines.
#' @export
plotPriorPost <- function(fittedModel, M=2e5, ci=.95){

  mfrow <- par()$mfrow
  hyp  <- hyp.eval <- fittedModel$mptInfo$hyperprior
  S <- length(fittedModel$mptInfo$thetaUnique)

  if(fittedModel$mptInfo$model == "traitMPT"){
    ww <- rWishart(n = M, df = fittedModel$mptInfo$hyperprior$df,
                   Sigma = fittedModel$mptInfo$hyperprior$V)
    ss <- array(apply(ww,3, solve), c(S,S,M))
    cc <- array(apply(ss,3, cov2cor), c(S,S,M))
  }

  for(s in 1:S){
    # both betaMPT and traitMPT have at least two hyperprior elements
    # that may have either length one (default) or length S:
    for(i in 1:2){
      if(is.na( hyp[[i]][s]))
        hyp[[i]][s] <- hyp[[i]][1]
      hyp.eval[[i]][s] <- sub("(", paste0("(",M,","),
                     sub("d", "r", hyp[[i]][s]),  fixed=TRUE)
    }
    aa <- eval(parse(text=hyp.eval[[1]][s]))
    bb <- eval(parse(text=hyp.eval[[2]][s]))

    if(fittedModel$mptInfo$model == "betaMPT"){
      # formulas for mean and SD of beta distribution:
      mean <- aa/(aa+bb)
      sd <-  sqrt(aa*bb/((aa+bb)^2*(aa+bb+1)))
      d.sd <- density(unlist(fittedModel$runjags$mcmc[,paste0("sd[",s,"]")]),
                      from=0, to=.5)
      ci.sd <- quantile(unlist(fittedModel$runjags$mcmc[,paste0("sd[",s,"]")]),
                        c((1-ci)/2,1-(1-ci)/2))

      prior.sd <- density(sd, from=0, to=.5)
      xlab.sd = "Group SD"
    }else{
      # probit transform:
      mean= sd <- pnorm(aa)
      sd <- bb * sqrt(ss[s,s,])

      d.sd <- density(unlist(fittedModel$runjags$mcmc[,paste0("sigma[",s,"]")]))
      prior.sd <- density(sd)
      xlab.sd = "Group SD (on latent probit scale)"
      ci.sd <- quantile(unlist(fittedModel$runjags$mcmc[,paste0("sigma[",s,"]")]),
                        c((1-ci)/2,1-(1-ci)/2))
    }
    d.mean <- density(unlist(fittedModel$runjags$mcmc[,paste0("mean[",s,"]")]), from=0, to=1)
    prior.mean <- density(mean, from=0, to=1)


    par(mfrow=1:2)
    tmp <- readline(prompt = "Press <Enter> to show the next plot.")
    plot(d.mean, main=paste0( "Group mean of ", fittedModel$mptInfo$thetaUnique[s]),
         xlab="Group mean")
    lines(prior.mean, col="blue", lty="dashed")
    abline(v= quantile(unlist(fittedModel$runjags$mcmc[,paste0("mean[",s,"]")]),
                               c((1-ci)/2,1-(1-ci)/2)), col="red", lty="dotted")
    plot(d.sd,   main=paste0("Group SD of ", fittedModel$mptInfo$thetaUnique[s]),
         xlab=xlab.sd)
    lines(prior.sd, col="blue", lty="dashed")
    abline(v=ci.sd,col="red", lty="dotted")
  }
  if(fittedModel$mptInfo$model == "traitMPT" & S>1){
    cnt <- 0
    for(s1 in 1:(S-1)){
      for(s2 in (s1+1):S){
        d.cor <- density(unlist(fittedModel$runjags$mcmc[,paste0("rho[",s1,",",s2,"]")]),
                         from=-1, to=1)
        prior.cor <- density(cc[s1,s2,], from=-1, to=1)
        if(cnt/2 == round(cnt/2))
          tmp <- readline(prompt = "Press <Enter> to show the next plot.")
        cnt <- cnt+1
        plot(d.cor, main=paste0( "Correlation between ",
                                  fittedModel$mptInfo$thetaUnique[s1], " and ",
                                  fittedModel$mptInfo$thetaUnique[s2]),
             xlab="Correlation (on latent probit scale)")
        lines(prior.cor, col="blue", lty="dashed")
        abline(v= quantile(unlist(fittedModel$runjags$mcmc[,paste0("rho[",s1,",",s2,"]")]),
                           c((1-ci)/2,1-(1-ci)/2)), col="red", lty="dotted")
      }
    }
  }

  par(mfrow=mfrow)

}

