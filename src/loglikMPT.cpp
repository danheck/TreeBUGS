// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/rmultinom.h>
// #include <RcppArmadilloExtensions/fixprob.h>
// #include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// MPT log-likelihood (without multinomial normalizing constants!)
// [[Rcpp::export]]
arma::vec loglikMPT(arma::mat theta,  // MPT parameters
                    arma::vec h,      // response frequencies
                    arma::mat a,
                    arma::mat b,
                    arma::vec c,
                    arma::vec map)
{
  // initialize predicted category probabilities
  arma::vec cat = arma::zeros<arma::vec>(h.n_elem);
  double branch;
  arma::vec ll(theta.n_rows);

  // loop across separate sets of parameters:
  for(arma::uword m=0; m<theta.n_rows; m++)
  {
    cat.zeros();
    // loop across MPT branches:
    for (arma::uword k=0;k<a.n_rows;k++)
    {
      branch = 1;
      // loop across MPT parameters
      for (arma::uword i=0; i<a.n_cols; i++)
      {
        branch *= pow(theta(m,i), a(k,i))* pow(1-theta(m,i), b(k,i));
      }
      branch *= c(k);           // constants
      cat(map(k)-1) += branch;  // mixture
    }
    ll(m) = arma::dot(log(cat), h);
    if(!arma::is_finite(ll(m)))
      ll(m) = - arma::datum::inf;
  }
  return ll;
}
