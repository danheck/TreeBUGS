
# prior:
#   (A) betaMPT: list(alpha="dgamma(1,.1)", beta="dgamma(1,.1)")
#   (B) traitMPT: list(mu="dnorm(0,1)", xi="dunif(0,10)", df=3, V=diag(2))
# M: number of samples
# S: optional, only for betaMPT: number of parameters
#' @importFrom MASS mvrnorm
sampleHyperprior <- function(prior, M, S){

  if(length(prior) == 4 && all(c("mu","xi","V","df") %in% names(prior))){
    prior <- prior[c("mu","xi","V","df")]
    model <- "traitMPT"
    S <- nrow(prior$V)
  }else if(length(prior) ==2 && all(c("alpha","beta") %in% names(prior))){
    prior <- prior[c("alpha","beta")]
    model <- "betaMPT"
    if(missing(S) || is.null(S)){
      S <- 1
    }
  }else{
    stop("Names of the list 'prior' do not match the hyperprior parameters of\n",
         "the traitMPT (xi, V, df) or betaMPT (alpha, beta)")
  }


  # both betaMPT and traitMPT have at least two hyperprior elements
  # that may have either length one (default) or length S:
  aa <- bb <- matrix(NA, M, S)
  hyp  <- hyp.eval <- prior
  for(s in 1:S){
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
    ## beta: alpha/beta  ; trait: mu/xi
    aa[,s] <- eval(parse(text=hyp.eval[[1]][s]))
    bb[,s] <- eval(parse(text=hyp.eval[[2]][s]))
  }

  if(model == "traitMPT"){
    ww <- rWishart(n = M, df = prior$df, Sigma = prior$V)
    ss <- array(apply(ww,3, solve), c(S,S,M))
    ###### SD + correlation:
    sig <- bb * sqrt(t(apply(ss, 3, diag)))
    rho <- array(apply(ss,3, cov2cor), c(S,S,M))
    ######### check with Var-Covar-Matrix:
    # Sigma <- ss
    # multi <- matrix(1, S, S)
    # for(mm in 1:M){
    #   diag(multi) <- bb[mm,]^2
    #   Sigma[,,mm] <- ss[,,mm] * multi
    # }
    # print(Sigma[,,50])
    # print((sig[50,] %*% t(sig[50,]) )* rho[,,50])
  }

  if(model == "betaMPT"){
    res <- list(alpha=aa, beta=bb)
  }else{
    res <- list(mu=aa, sigma=sig, rho=rho)# Sigma=Sigma)
  }
  attr(res, "model") <- model

  res
}

#
# debug(sampleHyperprior)
# sampleHyperprior(list(alpha="dgamma(1,.1)", beta="dgamma(1,.1)"), 100, S=3)
# sampleHyperprior(list(mu="dnorm(0,1)", xi="dunif(0,1)", df=3, V=diag(2)), M=50)
