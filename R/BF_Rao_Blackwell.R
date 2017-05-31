##########################################
### Barker & Link (2013):
### Rao-Blackwell estimate of posterior probabilities
#
#' Bayes Factors for Simple (Nonhierarchical) MPT Models
#'
#' Computes Bayes factors for simple (fixed-effects, nonhierarchical) MPT models with beta distributions as priors on the parameters.
#'
#' @param models list of models fitted with \code{\link{simpleMPT}}, e.g., \code{list(mod1, mod2)}
#' @param resample how much parameter posterior samples should be resampled per model
# @param scale scaling factor for the posterior approximation of the posterior samples
#' @inheritParams marginalMPT
#' @param store whether to save parameter samples
#' @details
#' Currently, this is only implemented for a single data set!
#'
#' Uses a Rao-Blackwellized version of the product-space method (Carlin & Chib, 1995) as proposed by Barker and Link (2013). First, posterior distributions of the MPT parameters are approximated by independent beta distributions. Second, for one a selected model, parameters are sampled from these proposal distributions. Third, the conditional probabilities to switch to a different model are computed and stored. Finally, the eigenvector with eigenvalue one of the matrix of switching probabilities provides an estimate of the posterior model probabilities.
#' @references
#' Barker, R. J., & Link, W. A. (2013). Bayesian multimodel inference by RJMCMC: A Gibbs sampling approach. The American Statistician, 67(3), 150-156.
#'
#' Carlin, B. P., & Chib, S. (1995). Bayesian model choice via Markov chain Monte Carlo methods. Journal of the Royal Statistical Society. Series B (Methodological), 57(3), 473-484.
#' @export
#' @seealso \code{\link{marginalMPT}}
BayesFactorMPT <- function(models,
                           #method="custom",
                           resample = 10000,
                           # scale=.1,
                           batches = 10,
                           store=FALSE,
                           nCPU = 1){

  if (!is.list(models) || any(sapply(models, class) != "simpleMPT"))
    stop("'models' must be a list models with fitted simpleMPT!")
  datas <- lapply(models, function(m) m$mptInfo$data)
  M <- length(models)
  if (M<2)
    stop("At least two models must be provided.")
  for (i in 2:M)
    if (any(datas[[1]] != datas[[i]]) )
      stop ("each model must have one vector of frequencies that must be identical for all ")
  if (nrow(datas[[1]]) != 1)
    warning("Only the first data set is used!")

  # 2. Approximate posteriors by beta densities
  betapars <- shape.prior <- list()
  for(m in 1:M){
    ab <- approximatePosterior(models[[m]], sample=2000)
    betapars[[m]] <- pmax(ab, 1) #pmax(ab*scale, 1)

    shape.prior[[m]] <- do.call("cbind", models[[m]]$mpt$hyperprior)
  }

  # 3. Loop 1 (rows of P): Model k
  #     => sample palette vector: phi = c(t1, ..., tk, ..., tm)
  #     => tk ~ resample from MCMC posteriors
  #     => ti ~ sample from approximation

  P <- array(NA, c(M,M,resample))

  # m = row of P (FROM which model to start jumping)
  row.P <- function(m){
    # sample proposal parameters
    theta <- vector("list", M)
    theta[[m]] <- resampling(models[[m]], resample=resample)
    theta[-m] <- lapply(betapars[-m],
                        function(bp) apply(bp, 1,
                                           function(s) rbeta(resample, s[1], s[2])))
    # loglik <- prior.pseudo <- prior.current <-
    posterior <- matrix(0, resample, M)

    # columns of P (TO which model to go)
    # for(m2 in 1:M){
    #   # sample proposal vector:
    #   if(m == m2){
    #     theta[[m]] <- resampling(models[[m]], resample=resample)
    #   }else{
    #     theta[[m2]] <- apply(betapars[[m2]], 1,
    #                          function(s) rbeta(resample, s[1], s[2]))
    #   }
    #   loglik[,m2] <- llMPT(theta[[m2]], mod = models[[m2]], id = 1)
    #   prior.pseudo[,m2] <- dProductBeta(x = theta[[m2]],
    #                                     shapes = betapars[[m2]])
    #   prior.current[,m2] <-
    #     dProductBeta(theta[[m2]], shapes = shape.prior[[m2]])
    # }

    loglik <- mapply(llMPT, pars = theta, mod = models, MoreArgs = list(id = 1))
    prior.pseudo <- mapply(dProductBeta, x = theta, shapes = betapars)
    prior.current <- mapply(dProductBeta, x = theta, shapes = shape.prior)

    ### 4. Loop 2 (entries in row k): Compute transition probabilities
    ###     => P(i|k) = P(y|tk,Mk) * P(tk|Mk) * prod(P(ti|Mk)) * P(Mk)
    # prior for full palette vector:
    posterior <-
      exp(loglik +                               # y | theta_k, M_k
            prior.current +                      # theta_k | M_k
            rowSums(prior.pseudo)-prior.pseudo)  # theta_i | M_k
    # for(m3 in 1:M){
    #   posterior[,m3] <-
    #     exp(loglik[,m3] +                                 # y | theta_k, M_k
    #           prior.current[,m3] +                        # theta_k | M_k
    #           rowSums(prior.pseudo[,-c(m3),drop=FALSE]))  # theta_i | M_k
    # }
    # P[m,,] <- t(posterior/rowSums(posterior))
    posterior/rowSums(posterior)
  }
  if(nCPU>1){
    cl <- makeCluster(nCPU)
    clusterExport(cl, c("resample", "M", "models", "betapars", "shape.prior"),
                  envir = environment())
    # tmp <- clusterEvalQ(cl, library(TreeBUGS))
    P.tmp <- parSapply(cl, 1:M, row.P, simplify = FALSE)
    stopCluster(cl)
  }else{
    P.tmp <- sapply(1:M, row.P, simplify = FALSE)
  }
  for(m in 1:M)
    P[m,,] <- t(P.tmp[[m]])
  rm(P.tmp) ; gc()

  #### 5. Rao-Blackwell Estimate: Average probabilities
  ####    Left eigenvector with eigenvalue = 1
  P.mean <- apply(P, 1:2, mean)
  ev2 <- Re(eigen(t(P.mean))$vec[,1])
  p.est <- ev2/sum(ev2)

  # batch estimate + SE
  idx <- rep(1:batches, each=round(resample/batches))
  if(length(idx) > resample) idx <- idx[1:resample]
  if(length(idx) < resample) idx <- c(idx, rep(batches,resample-length(idx)))
  tmp <- matrix(NA, batches, M)
  for(b in 1:batches){
    P.mean <- apply(P[,,idx == b], 1:2, mean)
    ev2 <- Re(eigen(t(P.mean))$vec[,1])
    tmp[b,] <- ev2/sum(ev2)
  }
  batchSE <- apply(tmp, 2, sd)/sqrt(batches)
  p.batch <- matrix(c(colMeans(tmp), batchSE), nrow=2, ncol=M,byrow = TRUE,
                    dimnames=list(c("Mean","SE"), Model=names(models)))

  dimnames(P.mean) <- list("from" = names(models),
                           "to" = names(models))
  names(p.est) <- names(models)
  res <- list("posterior" = p.est,
              "p.batch" = p.batch,
              "P.mean" = P.mean)
  if(store)
    res$samples <- P
  res
}


# compute priors for proposal probabilities (conditional on m3)
# sapply(theta, function(xx) xx[1,]) #### FIXED PROPOSAL!
# loglik[1,]    ###### conditionally independent
#
# #### conditonal:  psi | M_1
# pr1 <- (dProductBeta(theta[[1]], do.call("cbind", models[[1]]$mpt$hyperprior))[1]+
#           dProductBeta(x = theta[[2]], shapes = betapars[[2]])+
#           dProductBeta(x = theta[[3]], shapes = betapars[[3]]))[1]
# #### conditonal:  psi | M_2
# pr2 <- (dProductBeta(x = theta[[1]], shapes = betapars[[1]])+
#           dProductBeta(theta[[2]], do.call("cbind", models[[2]]$mpt$hyperprior))[1]+
#           dProductBeta(x = theta[[3]], shapes = betapars[[3]]))[1]
# #### conditonal:  psi | M_2
# pr3 <- (dProductBeta(x = theta[[1]], shapes = betapars[[1]])+
#           dProductBeta(x = theta[[2]], shapes = betapars[[2]]))[1]+
#   dProductBeta(theta[[3]], do.call("cbind", models[[3]]$mpt$hyperprior))[1]
# c(pr1, pr2, pr3) ########### prior conditional on each model
# posterior <- loglik[1] * c(pr1, pr2, pr3)
# posterior/sum(posterior)


# pp <- t(apply(P, 3, function(p)
#   Re(eigen(t(p))$vectors[,1])))
# pp <- pp/rowSums(pp)
