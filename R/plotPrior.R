
#' Plot Prior Distributions
#'
#' Plots prior distributions for group means, standard deviation, and correlations of MPT parameters across participants.
#' @param M number of random samples to approximate priors of group-level parameters
#' @param probitInverse which latent-probit parameters (for \code{\link{traitMPT}} model) to transform to probability scale. Either \code{"none"}, \code{"mean"} (simple transformation \eqn{\Phi(\mu)}), or \code{"mean_sd"} (see \code{\link{probitInverse}})
#' @inheritParams priorPredictive
#' @details This function samples from a set of hyperpriors (either for hierarchical traitMPT or betaMPT structure) to approximate the implied prior distributions on the parameters of interest (group-level mean, SD, and correlations of MPT parameters).
#' @export
#' @seealso \code{\link{priorPredictive}}
#' @examples
#' \dontrun{
#' # default priors for traitMPT:
#' plotPrior(list(mu = "dnorm(0,1)", xi="dunif(0,10)",
#'           V=diag(2), df=2+1), M=4000)
#'
#' # default priors for betaMPT:
#' plotPrior(list(alpha ="dgamma(1,.1)",
#'           beta = "dgamma(1,.1)"), M=4000)
#' }
# @importFrom evmix dbckden
plotPrior <- function(prior, probitInverse = "none", M=5000, nCPU=3){

  samples <- sampleHyperprior(prior, M, S=1, nCPU=nCPU)
  model <- attr(samples, "model")

  mfrow <- par()$mfrow
  qq <- seq(0,1, .05)
  bins <- min(60, round(M/40))
  histcol <- adjustcolor("gray", alpha.f = .7)

  if(model == "betaMPT"){
    S <- ncol(samples$alpha)
    # formulas for mean and SD of beta distribution:
    aa <- samples$alpha
    bb <- samples$beta
    mean <- aa/(aa+bb)
    sd <-  sqrt(aa*bb/((aa+bb)^2*(aa+bb+1)))
    # sd.sub <- sd[sample(M, min(5000,M))]
    # prior.sd <- density(sd, from=0, to=.5)
    # prior.sd <- dbckden(quantile(sd.sub, qq),  sd.sub,
    #                     lambda = .02, bcmethod = "beta2", xmax=.5)
  }else{
    # probit transform:
    S <- ncol(samples$mu)
    if(probitInverse %in% c("none","mean")){
      sd <- samples$sigma
      if(probitInverse == "none")
        mean <- samples$mu
      else
        mean <- pnorm(samples$mu)

      for(s in 1:S){
        sd[sd[,s]<=quantile(sd[,s], .997, na.rm = TRUE),s] <- NA ## remove extreme values
      }
    }
    if (probitInverse == "mean_sd"){
      mean <- sd <- samples$mu
      for(s in 1:S){
        mean_sd <- probitInverse(samples$mu[,s], samples$sigma[,s])
        mean[,s] <- mean_sd[,1]
        sd[,s] <- mean_sd[,2]
      }
    }

    # sd.sub <- sd[sample(M, min(5000,M))]
    # prior.sd <- density(sd, from = 0)
    # prior.sd <- dbckden(quantile(sd.sub, qq),  sd.sub, lambda = .2)#, bcmethod = "beta2", xmax=1)
  }
  # mean.sub <- mean[sample(M, min(5000,M))]
  # prior.mean <- dbckden(quantile(mean.sub, qq),  mean.sub,
  #                       lambda = .03, bcmethod = "beta2", xmax=1)


  par(mfrow=c(2,ifelse(model=="traitMPT",2,1)))
  hist(mean[,1], bins, freq=FALSE, col=histcol,
       main="Prior on group mean",
       xlab=paste0("Group mean",
                   ifelse(model=="traitMPT",
                          ifelse(probitInverse=="none"," (probit scale)", " (probability scale)"),
                          "")), border = histcol )
  # lines(quantile(mean.sub, qq), prior.mean)
  hist(sd[,1], bins, freq=FALSE, col=histcol, xlim=c(0, max(max(sd), .5)),
       main=paste0("Prior on group SD"),
       # ifelse(model=="traitMPT"," (on latent probit scale)","")),
       xlab=paste0("Group SD",
                   ifelse(model=="traitMPT",
                          ifelse(probitInverse=="mean_sd", " (probability scale)"," (probit scale)"),
                          "")), border = histcol )
  # lines(quantile(sd.sub, qq), prior.sd)
  if(S>1){
    if((model == "traitMPT" && max(length(prior$mu), length(prior$xi))>1) |
       (model == "betaMPT" && max(length(prior$alpha), length(prior$beta))>1  )){

      hist(mean[,s], bins, freq=FALSE, col=histcol,
           main=paste0("Prior on group mean ", s),
           xlab=paste0("Group mean" , ifelse(model=="traitMPT",
                              ifelse(probitInverse=="none"," (probit scale)", " (probability scale)"),
                              "")), border = histcol )
      hist(sd[,s], bins, freq=FALSE, col=histcol, xlim=c(0, max(max(sd), .5)),
           main=paste0("Prior on group SD ", s),
           xlab=paste0("Group SD",
                       ifelse(model=="traitMPT",
                              ifelse(probitInverse=="mean_sd", " (probability scale)"," (probit scale)"),
                              "")), border = histcol )
    }
  }

  if(model == "traitMPT" && S>1){
    for(s1 in 1:(S-1)){
      for(s2 in (s1+1):S){
        # rho.sub <- samples$rho[,,sample(M, min(5000,M))]
        # prior.cor <- dbckden(quantile(rho.sub[s1,s2,], qq)+1,
        #                      rho.sub[s1,s2,]+1,
        #                      lambda = .05, bcmethod = "beta2", xmax=2)
        hist(samples$rho[s1,s2,], bins, freq=FALSE, col=histcol,
             main=paste0("Correlation: ", s1, " and ",s2),
             xlab="Correlation (probit scale)", border = histcol )
        # lines(quantile(rho.sub[s1,s2,], qq), prior.cor)
      }
    }
  }

  par(mfrow=mfrow)

}



# if(length(prior) == 4 && all(c("mu","xi","V","df") %in% names(prior))){
#   prior <- prior[c("mu","xi","V","df")]
#   model <- "traitMPT"
#   S <- nrow(prior$V)
# }else if(length(prior) ==2 && all(c("alpha","beta") %in% names(prior))){
#   prior <- prior[c("alpha","beta")]
#   model <- "betaMPT"
#   S <- 2
# }else{
#   stop("Names of the list 'prior' do not match the hyperprior parameters of\n",
#        "the traitMPT (xi, V, df) or betaMPT (alpha, beta)")
# }

# if(model == "traitMPT"){
#   ww <- rWishart(n = M, df = prior$df, Sigma = prior$V)
#   ss <- array(apply(ww,3, solve), c(S,S,M))
#   cc <- array(apply(ss,3, cov2cor), c(S,S,M))
# }

# hyp  <- hyp.eval <- prior
# for(s in 1:1){
#   # both betaMPT and traitMPT have at least two hyperprior elements
#   # that may have either length one (default) or length S:
#   for(i in 1:2){
#     if(is.na( hyp[[i]][s]))
#       hyp[[i]][s] <- hyp[[i]][1]
#     hyp.eval[[i]][s] <- sub("(", paste0("(",M,","),
#                             sub("d", "r", hyp[[i]][s]),  fixed=TRUE)
#
#     # JAGS uses precision for normal distribution:
#     if(grepl("norm",hyp.eval[[i]][s] )){
#       tmp <- strsplit(hyp.eval[[i]][s], c( "[,\\(\\)]"), perl=F)[[1]]
#       tmp[4] <- 1/sqrt(as.numeric(tmp[4]))
#       hyp.eval[[i]][s] <- paste0(tmp[1],"(",tmp[2],",",tmp[3],",",tmp[4],")")
#     }
#   }
#   aa <- eval(parse(text=hyp.eval[[1]][s]))
#   bb <- eval(parse(text=hyp.eval[[2]][s]))
