% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BF_marginalMPT.R
\name{marginalMPT}
\alias{marginalMPT}
\title{Marginal Likelihood for Simple MPT}
\usage{
marginalMPT(
  eqnfile,
  data,
  restrictions,
  alpha = 1,
  beta = 1,
  dataset = 1,
  method = "importance",
  posterior = 500,
  mix = 0.05,
  scale = 0.9,
  samples = 10000,
  batches = 10,
  show = TRUE,
  cores = 1
)
}
\arguments{
\item{eqnfile}{The (relative or full) path to the file that specifies the MPT
model (standard .eqn syntax). Note that category labels must start with a
letter (different to multiTree) and match the column names of \code{data}.
Alternatively, the EQN-equations can be provided within R as a character
value (cf. \code{\link{readEQN}}). Note that the first line of an .eqn-file
is reserved for comments and always ignored.}

\item{data}{The (relative or full) path to the .csv file with the data (comma
separated; category labels in first row). Alternatively: a data frame or
matrix (rows=individuals, columns = individual category frequencies,
category labels as column names)}

\item{restrictions}{Specifies which parameters should be (a) constant (e.g.,
\code{"a=b=.5"}) or (b) constrained to be identical (e.g., \code{"Do=Dn"})
or (c) treated as fixed effects (i.e., identical for all participants;
\code{"a=b=FE"}). Either given as the path to a text file with restrictions
per row or as a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}.
Note that numbers in .eqn-equations (e.g., \code{d*(1-g)*.50}) are directly
interpreted as equality constraints.}

\item{alpha}{first shape parameter(s) for the beta prior-distribution of the
MPT parameters \eqn{\theta_s} (can be a named vector to use a different
prior for each MPT parameter)}

\item{beta}{second shape parameter(s)}

\item{dataset}{for which data set should Bayes factors be computed?}

\item{method}{either \code{"importance"} (importance sampling using a mixture
of uniform and beta-aproximation of the posterior) or \code{"prior"} (brute
force Monte Carlo sampling from prior)}

\item{posterior}{number of posterior samples used to approximate
importance-sampling densities (i.e., beta distributions)}

\item{mix}{mixture proportion of the uniform distribution for the
importance-sampling density}

\item{scale}{how much should posterior-beta approximations be downscaled to
get fatter importance-sampling density}

\item{samples}{total number of samples from parameter space}

\item{batches}{number of batches. Used to compute a standard error of the
estimate.}

\item{show}{whether to show progress}

\item{cores}{number of CPUs used}
}
\description{
Computes the marginal likelihood for simple (fixed-effects, nonhierarchical)
MPT models.
}
\details{
Currently, this is only implemented for a single data set!

If \code{method = "prior"}, a brute-force Monte Carlo method is used and
parameters are directly sampled from the prior.Then, the likelihood is
evaluated for these samples and averaged (fast, but inefficient).

Alternatively, an importance sampler is used if \code{method = "importance"},
and the posterior distributions of the MPT parameters are approximated by
independent beta distributions. Then each parameter \eqn{s} is sampled from
the importance density:

\eqn{mix*U(0,1) + (1-mix)*Beta(scale*a_s, scale*b_s)}
}
\examples{
# 2-High-Threshold Model
eqn <- "## 2HTM ##
   Target  Hit  d
   Target  Hit  (1-d)*g
   Target  Miss (1-d)*(1-g)
   Lure    FA   (1-d)*g
   Lure    CR   (1-d)*(1-g)
   Lure    CR   d"
data <- c(
  Hit = 46, Miss = 14,
  FA = 14, CR = 46
)

# weakly informative prior for guessing
aa <- c(d = 1, g = 2)
bb <- c(d = 1, g = 2)
curve(dbeta(x, aa["g"], bb["g"]))

# compute marginal likelihood
htm <- marginalMPT(eqn, data,
  alpha = aa, beta = bb,
  posterior = 200, samples = 1000
)
# second model: g=.50
htm.g50 <- marginalMPT(eqn, data, list("g=.5"),
  alpha = aa, beta = bb,
  posterior = 200, samples = 1000
)

# Bayes factor
# (per batch to get estimation error)
bf <- htm.g50$p.per.batch / htm$p.per.batch
mean(bf) # BF
sd(bf) / sqrt(length(bf)) # standard error of BF estimate

}
\references{
Vandekerckhove, J. S., Matzke, D., & Wagenmakers, E. (2015).
  Model comparison and the principle of parsimony. In Oxford Handbook of
  Computational and Mathematical Psychology (pp. 300-319). New York, NY:
  Oxford University Press.
}
\seealso{
\code{\link{BayesFactorMPT}}
}
