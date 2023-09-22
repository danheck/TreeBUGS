#' TreeBUGS: Hierarchical Multinomial Processing Tree Modeling
#'
#' @description{
#' \if{html}{\figure{TreeBUGS.png}{options: width='120' alt='logo' style='float: right'}}
#' \if{latex}{\figure{TreeBUGS.png}{options: width=0.5in}}
#'
#' Uses standard MPT files in the .eqn-format (Moshagen, 2010) to fit
#' hierarchical Bayesian MPT models. Note that the software JAGS is required
#' (\url{https://mcmc-jags.sourceforge.io/}). }
#'
#'
#' The core functions either fit a Beta-MPT model (\code{\link{betaMPT}}; Smith
#' & Batchelder, 2010) or a latent-trait MPT model (\code{\link{traitMPT}};
#' Klauer, 2010). A fitted model can be inspected using convenient summary and
#' plot functions tailored to hierarchical MPT models.
#'
#' Detailed explanations and examples can be found in the package vignette,
#' accessible via \code{vignette("TreeBUGS")}
#'
#' @author Daniel W. Heck, Denis Arnold, & Nina Arnold
#' @docType package
#'
#' @importFrom runjags run.jags extract autoextend.jags extend.jags
#' @importFrom coda gelman.diag effectiveSize as.mcmc.list as.mcmc
#' @importFrom utils read.csv write.table write.csv capture.output count.fields
#' @importFrom graphics axis plot points  segments abline boxplot curve hist
#'   lines par
#' @importFrom grDevices rainbow adjustcolor
#' @importFrom stats pnorm rnorm runif sd qnorm dnorm dbeta quantile rWishart
#'   ave pchisq window rbinom var
#' @importFrom parallel parSapply
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom MASS fitdistr
#' @useDynLib "TreeBUGS", .registration=TRUE
#'
#'
#' @section Citation:
#'
#' If you use TreeBUGS, please cite the software as follows:
#'
#' Heck, D. W., Arnold, N. R., & Arnold, D. (2018).
#' TreeBUGS: An R package for hierarchical multinomial-processing-tree modeling.
#' \emph{Behavior Research Methods, 50}, 264–284.
#' \doi{10.3758/s13428-017-0869-7}
#'
#'
#' @section Tutorial:
#'
#' For a tutorial on MPT modeling (including hierarchical modeling in TreeBUGS), see:
#'
#' Schmidt, O., Erdfelder, E., & Heck, D. W. (2023).
#' How to develop, test, and extend multinomial processing tree models: A tutorial.
#' \emph{Psychological Methods}.
#' \doi{10.1037/met0000561}.
#' (Preprint: \url{https://psyarxiv.com/gh8md/})
#'
#'
#' @references
#'
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models:
#' A latent-trait approach.
#' \emph{Psychometrika, 75}, 70-98.
#' \doi{10.1007/s11336-009-9141-0}
#'
#' Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015).
#' Bayesian estimation of multinomial processing tree models with heterogeneity
#' in participants and items.
#' \emph{Psychometrika, 80}, 205-235.
#' \doi{10.1007/s11336-013-9374-9}
#'
#' Moshagen, M. (2010).
#' multiTree: A computer program for the analysis of multinomial processing
#' tree models.
#' \emph{Behavior Research Methods, 42}, 42-54.
#' \doi{10.3758/BRM.42.1.42}
#'
#' Smith, J. B., & Batchelder, W. H. (2008).
#' Assessing individual differences in categorical data.
#' \emph{Psychonomic Bulletin & Review, 15}, 713-731.
#' \doi{10.3758/PBR.15.4.713}
#'
#' Smith, J. B., & Batchelder, W. H. (2010).
#' Beta-MPT: Multinomial processing tree models for addressing
#' individual differences.
#' \emph{Journal of Mathematical Psychology, 54}, 167-183.
#' \doi{10.1016/j.jmp.2009.06.007}
#'
"_PACKAGE"
