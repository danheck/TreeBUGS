
#' Plot Prior Distributions
#'
#' Plots prior distributions for group means, standard deviation, and correlations of MPT parameters across participants.
#'
#' @param prior a named list defining the priors. For the \link{traitMPT}, the default is \code{list(mu = "dnorm(0,1)", xi="dunif(0,10)", V=diag(S), df=S+1)}, where S is the number of free parameters. For the \link{betaMPT}, the default is \code{list(alpha ="dgamma(1,.1)", beta = "dgamma(1,.1)")}
# @param S number of free parameters (only relevant for inverse Wishart prior of latent-trait MPT)
#' @param M number of random samples to approximate prior distributions
#' @details This function samples from a set of hyperpriors (either for hierarchical traitMPT or betaMPT structure) to approximate the implied prior distributions on the parameters of interest (mean, SD, correlations of MPT parameters).
#' @export
#' @examples
#' # default priors for traitMPT:
#' plotPrior(list(mu = "dnorm(0,1)", xi="dunif(0,10)",
#'                V=diag(2), df=2+1), M=4000)
#'
#' # default priors for betaMPT:
#' plotPrior( list(alpha ="dgamma(1,.1)",
#'                 beta = "dgamma(1,.1)"), M=4000)
#' @importFrom evmix dbckden
plotPrior <- function(prior, M=10000){

  mfrow <- par()$mfrow
  if(length(prior) == 4 && all(c("mu","xi","V","df") %in% names(prior))){
    prior <- prior[c("mu","xi","V","df")]
    model <- "traitMPT"
    S <- nrow(prior$V)
  }else if(length(prior) ==2 && all(c("alpha","beta") %in% names(prior))){
    prior <- prior[c("alpha","beta")]
    model <- "betaMPT"
    S <- 2
  }else{
    stop("Names of the list 'prior' do not match the hyperprior parameters of\n",
         "the traitMPT (xi, V, df) or betaMPT (alpha, beta)")
  }

  if(model == "traitMPT"){
    ww <- rWishart(n = M, df = prior$df, Sigma = prior$V)
    ss <- array(apply(ww,3, solve), c(S,S,M))
    cc <- array(apply(ss,3, cov2cor), c(S,S,M))
  }

  qq <- seq(0,1, .05)
  bins <- min(80, round(M/40))
  histcol <- adjustcolor("gray", alpha.f = .7)
  hyp  <- hyp.eval <- prior
  for(s in 1:1){
    # both betaMPT and traitMPT have at least two hyperprior elements
    # that may have either length one (default) or length S:
    for(i in 1:2){
      if(is.na( hyp[[i]][s]))
        hyp[[i]][s] <- hyp[[i]][1]
      hyp.eval[[i]][s] <- sub("(", paste0("(",M,","),
                              sub("d", "r", hyp[[i]][s]),  fixed=TRUE)

      # JAGS uses precision for normal distribution:
      if(grepl("norm",hyp.eval[[i]][s] )){
        tmp <- strsplit(hyp.eval[[i]][s], c( "[,\\(\\)]"), perl=F)[[1]]
        tmp[4] <- 1/sqrt(as.numeric(tmp[4]))
        hyp.eval[[i]][s] <- paste0(tmp[1],"(",tmp[2],",",tmp[3],",",tmp[4],")")
      }
    }
    aa <- eval(parse(text=hyp.eval[[1]][s]))
    bb <- eval(parse(text=hyp.eval[[2]][s]))

    if(model == "betaMPT"){
      # formulas for mean and SD of beta distribution:
      mean <- aa/(aa+bb)
      sd <-  sqrt(aa*bb/((aa+bb)^2*(aa+bb+1)))
      prior.sd <- density(sd, from=0, to=.5)
      xlab.sd = "Group SD"
    }else{
      # probit transform:
      mean <- pnorm(aa)
      sd <- bb * sqrt(ss[s,s,])   # scaled SD on group level
      sd <- sd[sd<=quantile(sd, .997)]
      prior.sd <- density(sd, from = 0)
      xlab.sd = "Group SD (on latent probit scale)"
    }
    # prior.mean <- density(mean, from=0, to=1, kernel="cosine")
    prior.mean <- dbckden(quantile(mean, qq),  mean, lambda = .03, bcmethod = "beta2", xmax=1)


    par(mfrow=c(2,2))
    # tmp <- readline(prompt = "Press <Enter> to show the next plot.")
    hist(mean, bins, freq=FALSE, col=histcol,
         main=" Prior on group mean", xlab="Group mean", border = histcol )
    lines(quantile(mean, qq), prior.mean)
    hist(sd, bins, freq=FALSE, col=histcol, xlim=c(0, max(quantile(sd, .995), .5)),
         main=paste0("Prior on SD",
                     ifelse(model=="traitMPT"," (on latent probit-scale)","")),
         xlab="SD", border = histcol )
    lines(prior.sd)
  }
  if(model == "traitMPT" & S>1){
    cnt <- 0
    for(s1 in 1:(S-1)){
      for(s2 in (s1+1):S){
        prior.cor <- dbckden(quantile(cc[s1,s2,], qq)+1,  cc[s1,s2,]+1,
                              lambda = .03, bcmethod = "beta2", xmax=2)
        # prior.cor <- density(cc[s1,s2,], from=-1, to=1)
        if(cnt/2 == round(cnt/2))
          # tmp <- readline(prompt = "Press <Enter> to show the next plot.")
        cnt <- cnt+1
        hist(cc[s1,s2,], bins, freq=FALSE, col=histcol,
             main=paste0("Correlation: ", s1, " and ",s2),
             xlab="Correlation (on latent probit scale)", border = histcol )
        lines(quantile(cc[s1,s2,], qq), prior.cor)
      }
    }
  }

  par(mfrow=mfrow)

}

