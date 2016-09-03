#' TreeBUGS: Hierarchical Multinomial Processing Tree Modeling
#'
#' Uses standard MPT files in the .eqn-format (Moshagen, 2010) to fit hierarchical Bayesian MPT models. Note that the software JAGS is required (\url{http://mcmc-jags.sourceforge.net}).
#'
#' The core functions either fit a Beta-MPT model (\code{\link{betaMPT}};Smith & Batchelder, 2010) or a latent-trait MPT model (\code{\link{traitMPT}}; Klauer, 2010). A fitted model can be inspected using convenient summary and plot functions tailored to hierarchical MPT models.
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("TreeBUGS")}
#' @author Daniel Heck, Denis Arnold, & Nina Arnold
#' @docType package
#' @name TreeBUGS
#' @importFrom runjags run.jags extract autoextend.jags extend.jags
#run.jags  summary.runjags extract autorun.jags autoextend.jags
#' @importFrom coda gelman.diag effectiveSize as.mcmc.list as.mcmc
#' @importFrom utils read.csv write.table write.csv
#' @importFrom graphics axis plot points  segments abline boxplot curve hist lines par
#' @importFrom grDevices rainbow adjustcolor
#' @importFrom stats pnorm rnorm runif sd qnorm dnorm dbeta quantile rWishart ave pchisq
#' @importFrom parallel parSapply
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib TreeBUGS
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#'
#' Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015). Bayesian estimation of multinomial processing tree models with heterogeneity in participants and items. Psychometrika, 80, 205-235.
#'
#' Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. Behavior Research Methods, 42, 42-54.
#'
#'
#' Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167-183.
#'
#'
"_PACKAGE"


