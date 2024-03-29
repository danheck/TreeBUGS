% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WAIC.R
\name{WAIC}
\alias{WAIC}
\alias{print.waic}
\alias{print.waic_difference}
\alias{-.waic}
\title{WAIC: Widely Applicable Information Criterion}
\usage{
WAIC(
  fittedModel,
  n.adapt = 1000,
  n.chains = 3,
  n.iter = 10000,
  n.thin = 1,
  summarize = FALSE
)

\method{print}{waic}(x, ...)

\method{print}{waic_difference}(x, ...)

\method{-}{waic}(e1, e2)
}
\arguments{
\item{fittedModel}{fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})}

\item{n.adapt}{number of adaptation samples.}

\item{n.chains}{number of chains (no parallel computation).}

\item{n.iter}{number of iterations after burnin.}

\item{n.thin}{Thinning rate.}

\item{summarize}{deprecated argument only available for backwards compatibility}

\item{x}{An object of class \code{waic} or \code{waic_difference} to be printed.}

\item{...}{Further arguments that may be passed to print methods.}

\item{e1, e2}{Two objects of class \code{waic} to be compared.}
}
\value{
Function \code{WAIC()} returns an object of class \code{waic}, which is basically
a list containing three vectors \code{p_waic}, \code{deviance}, and \code{waic}, with
separate values for each observed node
(i.e., for all combinations of persons and free categories).

For these objects, a \code{print()} method exists, which
also calculates the standard error of the estimate of WAIC.

For backwards compatibility, if \code{WAIC()} is called with \code{summarize = TRUE},
a vector with values \code{p_waic}, \code{deviance}, \code{waic}, and \code{se_waic} is returned.

WAIC values from two models can be compared by using the \code{-} operator;
the result is an object of class \code{waic_difference}.
}
\description{
Implementation of the WAIC for model comparison.
}
\details{
WAIC provides an approximation of predictive accuracy with respect
to out-of-sample deviance. The uncertainty of the WAIC for the given number
of observed nodes (i.e., number of free categories times the number of
participants) is quantified by the standard error of WAIC \code{"se_waic"}
(cf. Vehtari et al., 2017). In contrast, to assess whether the approximation
uncertainty due to MCMC sampling (not sample size) is sufficiently low, it is
a good idea to fit each model twice and compute WAIC again to assess the
stability of the WAIC values.

For more details, see Vehtari et al. (2017) and the following discussion
about the JAGS implementation (which is currently an experimental feature of
JAGS 4.3.0):

\url{https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/}
}
\examples{
\dontrun{

#### WAIC for a latent-trait MPT model:
fit <- traitMPT(...)
WAIC(fit)


#### pairwise comparison of two models:

# (1) compute WAIC per model
waic1 <- WAIC(fit1)
waic2 <- WAIC(fit2)

# (2) WAIC difference
waic1 - waic2
}


}
\references{
Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian
model evaluation using leave-one-out cross-validation and WAIC. Statistics
and Computing, 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
}
