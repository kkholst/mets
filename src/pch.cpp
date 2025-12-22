// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo/Lighter>
#include <cmath>

using namespace Rcpp;
using arma::datum;

// [[Rcpp::export(name = ".rpch")]]
arma::vec rpch(unsigned n, std::vector<double> lambda, std::vector<double> time) {
  auto K = lambda.size();
  arma::vec res = -log(runif(n))/lambda[0] + time[0];
  for (unsigned i=0; i<n; i++) {
    for (unsigned k=1; k<K; k++) {
      if (res(i)<time[k]) break;
	     else res(i) = -log(R::runif(0,1))/lambda[k]+time[k]; // unif_rand()
    }
  }
  return ( res );
}

// [[Rcpp::export(name = ".cpch")]]
arma::vec cpch(arma::vec &x, std::vector<double> lambda, std::vector<double> time) {
  auto K = lambda.size();
  auto n = x.size();
  arma::vec res(n); res.fill(0);
  for (unsigned k=0; k<K; k++) {
    arma::uvec ind = (x >= time[k]);
    for (unsigned i=0; i<n; i++) {
      if (ind[i]>0) {
	double t0 = std::fmin(x[i]-time[k], time[k+1]-time[k]);
	res[i] += lambda[k]*t0;
      }
    }
  }
  return ( res );
}
