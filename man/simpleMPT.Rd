% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSimpleMPT.R
\name{simpleMPT}
\alias{simpleMPT}
\title{C++ Sampler for Standard (Nonhierarchical) MPT Models}
\usage{
simpleMPT(
  eqnfile,
  data,
  restrictions,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 3,
  n.chains = 3,
  ppp = 0,
  alpha = 1,
  beta = 1,
  parEstFile,
  posteriorFile,
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

\item{n.iter}{Number of iterations per chain (including burnin samples). See
\code{\link[runjags]{run.jags}} for details.}

\item{n.burnin}{Number of samples for burnin (samples will not be stored and
removed from n.iter)}

\item{n.thin}{Thinning rate.}

\item{n.chains}{number of MCMC chains (sampled in parallel, which can be
changed via the additional argument \code{n.sim = 1}).}

\item{ppp}{number of samples to compute  posterior predictive p-value (see
\code{\link{posteriorPredictive}})}

\item{alpha}{first shape parameter(s) for the beta prior-distribution of the
MPT parameters \eqn{\theta_s} (can be a named vector to use a different
prior for each MPT parameter)}

\item{beta}{second shape parameter(s)}

\item{parEstFile}{Name of the file to with the estimates should be stored
(e.g., "parEstFile.txt")}

\item{posteriorFile}{path to RData-file where to save the model including
MCMC posterior samples (an object named \code{fittedModel}; e.g.,
\code{posteriorFile="mcmc.RData"})}

\item{cores}{number of CPUs to be used}
}
\description{
Fast Gibbs sampler in C++ that is tailored to the standard fixed-effects MPT
model (i.e., fixed-effects, non-hierarchical MPT). Assumes independent
parameters per person if a matrix of frequencies per person is supplied.
}
\details{
Beta distributions with fixed shape parameters \eqn{\alpha} and
  \eqn{\beta} are used. The default \eqn{\alpha=1} and \eqn{\beta=1} assumes
  uniform priors for all MPT parameters.
}
\examples{
\dontrun{
# fit nonhierarchical MPT model for aggregated data (see ?arnold2013):
EQNfile <- system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS")
d.encoding <- subset(arnold2013, group == "encoding", select = -(1:4))
fit <- simpleMPT(EQNfile, colSums(d.encoding),
  restrictions = list("D1=D2=D3", "d1=d2", "a=g")
)
# convergence
plot(fit)
summary(fit)
}
}
\author{
Daniel Heck
}
