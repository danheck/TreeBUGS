// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/rmultinom.h>
#include <RcppArmadilloExtensions/fixprob.h>

using namespace Rcpp;


/* ############################ */
// CONDITIONAL POSTERIOR for alpha-parameter

// (for beta: use (1-theta) !!)
// x: previous value of parameter (e.g., to sample x=alpha[m-1])
// fixed: other parameter (e.g., conditional on fixed=beta[m])
// theta: vector of individual theta-parameters
// shape, rate: hyperprior for alpha, beta ~ Gamma(shape, rate)
double postAB(double x,
              double fixed,
              NumericVector theta,
              double shape=1.0,
              double rate=0.1){

  return R::dgamma(x, shape, 1/(rate-sum(log(theta))), 1)- theta.size() * R::lbeta( x, fixed) ;
}


/* ############################ */
// SLICE SAMPLER for alpha (for beta: use (1-theta)!):

// x: previous value of parameter (e.g., to sample x=alpha[m-1])
// fixed: other parameter (e.g., conditional on fixed=beta[m])
// theta: vector of individual theta-parameters
// shape, rate: hyperprior for alpha, beta ~ Gamma(shape, rate)
// eps, steps: precision/maximum number of steps to find upper/lower bound of slice
// **[[Rcpp::export]]
double sliceAB(double x,
               double fixed,
               NumericVector theta,
               double shape=1.0,
               double rate=0.1,
               double eps =0.01,
               int steps=5000){

  // get height logfx at x-coordinate:
  double logfx = postAB(x, fixed, theta, shape, rate);
  double z = logfx - R::rexp(1) ;

  // get LOWER x-axis boundary for intersection with z (=> bisection method)
  arma::vec xx,yy = arma::zeros(3);
  xx << 0.001 << x/2 << x;
  for(int i=0;i<3;i++){
    yy(i) = postAB(xx(i), fixed, theta, shape, rate) - z;
  }

  int cnt = 0;
  while( (norm(yy) > eps)  & (cnt < steps) ){
    cnt++;
    if( ( (yy(0)<0) & (yy(1)>0) ) | ( (yy(0)>0) & (yy(1)<0) )){
      xx(2) = xx(1);
    }else if( ( (yy(2)<0) & (yy(1)>0) ) | ( (yy(2)>0) & (yy(1)<0) )){
      xx(0) = xx(1);
    }

    // bisection proposal (midpoint of interval):
    xx(1) = (xx(2)-xx(0))/2 + xx(0);
    // secant method: might result in negative proposals!

    for(int i=0;i<3;i++){
      yy(i) = postAB(xx(i), fixed, theta, shape, rate) - z;
    }
  }
  double lower = xx(1);

  // get UPPER x-axis boundary for intersection with z (=> bisection method)
  xx << x << x+(max(NumericVector::create(3000, x*5))+x*9)/10 <<
    max(NumericVector::create(3000, x*5));
  for(int i=0;i<3;i++){
    yy(i) = postAB(xx(i), fixed, theta, shape, rate) - z;
  }

  cnt = 0;
  while( (norm(yy) > eps)  & (cnt < steps) ){
    cnt++;
    if( ( (yy(0)<0) & (yy(1)>0) ) | ( (yy(0)>0) & (yy(1)<0) )){
      xx(2) = xx(1);
    }else if( ( (yy(2)<0) & (yy(1)>0) ) | ( (yy(2)>0) & (yy(1)<0) )){
      xx(0) = xx(1);
    }
    xx(1) = (xx(2)-xx(0))/2 + xx(0);
    for(int i=0;i<3;i++){
      yy(i) = postAB(xx(i), fixed, theta, shape, rate) - z;
    }
  }
  double upper = xx(1);

  // Rcout << "lower/upper" << lower << " " << upper;
  double sample = R::runif(lower, upper);
  return sample;
}



// [[Rcpp::export]]
List betampt(int M,
             int L,
             int nthin,
             arma::mat H,
             arma::mat a,
             arma::mat b,
             arma::vec c,
             arma::vec map,
             arma::vec shape,
             arma::vec rate){

  int S = a.n_cols;     // parameters
  int N = H.n_rows;     // participants
  int B = a.n_rows;     // branches
  arma::vec cats = unique(map);
  int K = cats.n_elem;  // categories

  // MCMC matrices:
  arma::mat alpha(M, S);
  arma::mat beta(M, S);

  arma::mat mu(M, S);
  arma::mat sig(M, S);
  arma::cube theta(M, N, S);

  // temporary MCMC objects:
  arma::mat Hfull =  arma::zeros<arma::mat>(N, B);
  arma::uvec sel;
  arma::vec fixprob;
  NumericVector fp;
  IntegerVector samp;
  NumericVector tvec;

  arma::mat p(N,S);
  arma::mat q(N,S);
  arma::vec br(B);

  arma::mat theta_m(N, S);
  arma::vec alpha_m(S);
  arma::vec beta_m(S);

  // initialize temporary MCMC objects with random numbers from [0, 1]:
  theta_m.randu();
  alpha_m.randu();
  beta_m.randu();


  //  initialize full-data matrix Hfull
  for(int n=0; n<N; n++){
    for(int k=0; k<K; k++){

      sel = find(map == k+1);
      fixprob = arma::ones<arma::vec>(sel.n_elem) /sel.n_elem ;
      fp = Rcpp::as<Rcpp::NumericVector>(wrap(fixprob));

      if(fixprob.n_elem == 1){
        Hfull(n,as_scalar(sel)) = H(n,k);
      }else{
        samp = RcppArmadillo::rmultinom(H(n,k), fp);
        for(arma::uword tt=0; tt<fixprob.n_elem; tt++){
          Hfull(n,arma::as_scalar(sel(tt))) = samp[tt];
        }
      }
      // Rcout << "Hfull=\n" << Hfull.row(n).t() << "H =\n" << H.row(n).t() << "\n\n";
    }
  }

  // MCMC loop ----
  for(int m=-L; m<M; m++){

    // thinning loop ----
    for (int thin=0; thin<nthin; thin++) {

      //  hierarchical part ----

      for(int s=0; s<S; s++){
        tvec = as<Rcpp::NumericVector>(wrap(reshape(theta_m.col(s),N,1)));
        alpha_m(s) = sliceAB(alpha_m(s), beta_m(s), tvec, shape(s), rate(s), .001);
        beta_m(s) = sliceAB(beta_m(s), alpha_m(s), 1-tvec, shape(s), rate(s), .001);
      }


      // MPT part ----

      // # 3. Sample θ from θ|σ, µ, Hfull
      for(int s=0; s<S; s++){
        for(int n=0; n<N; n++){
          p(n,s) = dot(a.col(s), Hfull.row(n) );
          q(n,s) = dot(b.col(s), Hfull.row(n) );
          theta_m(n,s) = R::rbeta(p(n,s) + alpha_m(s), q(n,s) + beta_m(s)) ;
        }
      }

      // # 4. Sample Hfull from Hfull |θ, H
      for(int n=0; n<N; n++){

        br = c;
        for(int bb=0; bb<B; bb++){
          for(int s=0; s<S; s++){
            br(bb) = br(bb)* pow(theta_m(n,s), a(bb,s))*pow(1-theta_m(n,s), b(bb,s)) ;
          }
        }
        // Rcout << "theta=" << theta(span(m,m),span(n,n),span(0,S-1)) << "\nbr=" << br;
        for(int k=0; k<K; k++){

          sel = find(map == k+1);
          fixprob = br.elem(sel) /sum(br.elem(sel)) ;
          fp = as<Rcpp::NumericVector>(wrap(fixprob));
          // Rcout << "k=" << k << " fixprob=\n" << fixprob.t() << "Hfull=\n" << Hfull.row(n).t() << "\n\n";

          if(sel.n_elem == 1){
            Hfull(n,as_scalar(sel)) = H(n,k);
          }else{
            samp = RcppArmadillo::rmultinom(H(n,k), fp);
            for(arma::uword tt=0; tt<sel.n_elem; tt++){
              Hfull(n,arma::as_scalar(sel(tt))) = samp[tt];
            }
          }
        }

      }
      // Rcout << "Hfull=\n" << Hfull.row(n).t() << "H =\n" << H.row(n).t() << "\n\n";
      // Rcout << "\Here, Hfull=\n" << Hfull;
    } // end of thinning loop

    // Generate quantities only for stored iterations ----
    if(m >= 0) {
      theta.row(m) = theta_m;
      alpha.row(m) = alpha_m.t();
      beta.row(m) = beta_m.t();

      mu.row(m) = alpha.row(m)/(alpha.row(m) + beta.row(m));
      sig.row(m) = sqrt(mu.row(m) % (1 - mu.row(m)) / (alpha.row(m) + beta.row(m)+1));
    }
  } // end of MCMC loop

  return Rcpp::List::create(Rcpp::Named("mean") = mu,
                            Rcpp::Named("sd") = sig,
                            Rcpp::Named("alph") = alpha,
                            Rcpp::Named("bet") = beta,
                            Rcpp::Named("theta") = theta);
}





// [[Rcpp::export]]
List simplempt(int M,
               int L,
               int nthin,
               arma::mat H,
               arma::mat a,
               arma::mat b,
               arma::vec c,
               arma::vec map,
               arma::vec alpha,
               arma::vec beta){

  int S = a.n_cols;     // parameters
  int N = H.n_rows;     // participants
  int B = a.n_rows;     // branches
  arma::vec cats = unique(map);
  int K = cats.n_elem;  // categories

  //  Initialize MCMC array
  arma::cube theta(M, N, S);
  arma::mat theta_m(N, S);
  theta_m.randu();

  // temporary MCMC objects:
  arma::mat Hfull =  arma::zeros<arma::mat>(N, B);
  arma::uvec sel;
  arma::vec fixprob;
  NumericVector fp;
  IntegerVector samp;

  arma::mat p(N,S);
  arma::mat q(N,S);
  arma::vec br(B);

  //  initialize full-data matrix Hfull
  for(int n=0; n<N; n++){
    for(int k=0; k<K; k++){

      sel = find(map == k+1);
      fixprob = arma::ones<arma::vec>(sel.n_elem) /sel.n_elem ;
      fp = Rcpp::as<Rcpp::NumericVector>(wrap(fixprob));

      if(fixprob.n_elem == 1){
        Hfull(n,as_scalar(sel)) = H(n,k);
      }else{
        samp = RcppArmadillo::rmultinom(H(n,k), fp);
        for(arma::uword tt=0; tt<fixprob.n_elem; tt++){
          Hfull(n,arma::as_scalar(sel(tt))) = samp[tt];
        }
      }
      // Rcout << "Hfull=\n" << Hfull.row(n).t() << "H =\n" << H.row(n).t() << "\n\n";
    }
  }

  // ################################ MCMC loop
  for(int m=-L; m<M; m++){ // negative index for burnin, non-negative index for sampling
    for(int thin=0; thin<nthin; thin++) {

      // ################################ MPT part
      // # 3. Sample θ from θ|σ, µ, Hfull
      for(int s=0; s<S; s++){
        for(int n=0; n<N; n++){
          p(n,s) = dot(a.col(s), Hfull.row(n) );
          q(n,s) = dot(b.col(s), Hfull.row(n) );
          theta_m(n,s) = R::rbeta(p(n,s)+alpha(s), q(n,s) + beta(s)) ;
        }
      }

      // # 4. Sample Hfull from Hfull |θ, H
      for(int n=0; n<N; n++){

        br = c;
        for(int bb=0; bb<B; bb++){
          for(int s=0; s<S; s++){
            br(bb) = br(bb)* pow(theta_m(n,s), a(bb,s))*pow(1-theta_m(n,s), b(bb,s)) ;
          }
        }
        // Rcout << "theta=" << theta(span(m,m),span(n,n),span(0,S-1)) << "\nbr=" << br;
        for(int k=0; k<K; k++){

          sel = find(map == k+1);
          fixprob = br.elem(sel) /sum(br.elem(sel)) ;
          fp = as<Rcpp::NumericVector>(wrap(fixprob));
          // Rcout << "k=" << k << " fixprob=\n" << fixprob.t() << "Hfull=\n" << Hfull.row(n).t() << "\n\n";

          if(sel.n_elem == 1){
            Hfull(n,as_scalar(sel)) = H(n,k);
          }else{
            samp = RcppArmadillo::rmultinom(H(n,k), fp);  // Rcpp::
            for(arma::uword tt=0; tt<sel.n_elem; tt++){
              Hfull(n,arma::as_scalar(sel(tt))) = samp[tt];
            }
          }

        }
        // Rcout << "Hfull=\n" << Hfull.row(n).t() << "H =\n" << H.row(n).t() << "\n\n";
        // Rcout << "\Here, Hfull=\n" << Hfull;
      }
    }
    if(m >= 0) {
      theta.row(m) = theta_m;
    }
  }

  return Rcpp::List::create(Rcpp::Named("theta") = theta,
                            Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta);
}
