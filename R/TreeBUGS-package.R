#' TreeBUGS: Hierarchical MPT Modeling using BUGS
#'
#' Uses standard MPT files in the .eqn-format (Moshagen, 2010) to fit hierarchical Bayesian MPT models. Note that the software JAGS is required (\url{http://mcmc-jags.sourceforge.net}).
#'
#' The core functions either fit a Beta-MPT model (\code{\link{betaMPT}};Smith & Batchelder, 2010) or a latent-trait MPT model (\code{\link{traitMPT}}; Klauer, 2010). A fitted model can be inspected using convenient summary and plot functions tailored to hierarchical MPT models.
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("TreeBUGS")}
#' @author Daniel Heck, Denis Arnold, & Nina Arnold
#' @docType package
#' @name TreeBUGS
#' @import R2WinBUGS
#' @import R2jags
#' @importFrom utils read.csv write.table
#' @importFrom graphics axis plot points  segments
#' @importFrom grDevices rainbow
#' @importFrom stats pnorm rnorm runif sd qnorm
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


