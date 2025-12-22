/*!
  @file utils_interface.cpp
  @author Klaus K. Holst, Thomas Scheike
  @copyright 2025

  @brief Weibull likelihood and score

  The relevant bindings are created in \c RcppExports.cpp, \c RcppExports.h
  with \c Rcpp::compileAttributes()
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo/Lighter>

arma::vec chazw(const arma::vec &time, const arma::vec &rate,
                const arma::vec &shape) {
    return rate % arma::pow(time, shape);
}

using namespace Rcpp;

// [[Rcpp::export(name=".logl_weibull")]]
double logl_weibull(const arma::vec &par, // parameter vector
                    const arma::vec &entry,   // entry time
                    const arma::vec &exit,    // exit time
                    const arma::vec &status, // censoring (false is censored)
                    const arma::mat &X, // design matrix for log-rate
                    const arma::mat &Z // design matrix for log-shape
) {

    unsigned np1 = X.n_cols;
    unsigned np2 = Z.n_cols;
    arma::colvec par1(np1);
    arma::colvec par2(np2);;
    for (unsigned i = 0; i < np1; i++) par1[i] = par[i];
    for (unsigned i = 0; i < np2; i++) par2[i] = par[i + np1];
    arma::vec lp = X * par1;
    arma::vec rate = arma::exp(lp);
    arma::vec shape = arma::exp(Z * par2);
    arma::vec loghaz = lp + arma::log(shape) +
        (shape - 1) % arma::log(exit);
    arma::vec logl = status % loghaz -
        chazw(exit, rate, shape) + chazw(entry, rate, shape);
    return sum(logl);
}


arma::vec dchazw1(const arma::vec &time, const arma::vec &rate,
                  const arma::vec &shape) {
    return arma::pow(time, shape);
}


// [[Rcpp::export(name=".score_weibull")]]
arma::mat score_weibull(const arma::vec &par,    // parameter vector
              const arma::vec &entry,  // entry time
              const arma::vec &exit,   // exit time
              const arma::vec &status, // censoring (false is censored)
              const arma::mat &X,      // design matrix for log-rate
              const arma::mat &Z,      // design matrix for log-shape
              bool indiv = true
) {

    unsigned np1 = X.n_cols;
    unsigned np2 = Z.n_cols;
    arma::colvec par1(np1);
    arma::colvec par2(np2);;
    for (unsigned i = 0; i < np1; i++) par1[i] = par[i];
    for (unsigned i = 0; i < np2; i++) par2[i] = par[i + np1];
    arma::vec lp = X * par1;
    arma::vec rate = arma::exp(lp);
    arma::vec shape = arma::exp(Z * par2);

    arma::vec dloghaz1 = status / rate; // derivatitve of log-hazard wrt rate
    arma::vec dloghaz2 = status / shape + status % arma::log(exit); // derivatitve wrt shape
    arma::vec dchaz1entry = dchazw1(entry, rate, shape); // same for cum.haz evaluated in entry & exit time
    arma::vec dchaz1exit = dchazw1(exit, rate, shape);
    arma::vec logentry(entry.n_elem);
    arma::vec logexit(exit.n_elem);
    // res[idx0] < -0
    for (unsigned i = 0; i < entry.n_elem; i++) {
        if (entry[i] == 0) {
            logentry[i] = 0; //.VolumeIcon.icns conv. log(x)x = 0 when x = 0
        } else logentry[i] = std::log(entry[i]);
        if (exit[i] == 0) {
            logexit[i] = 0;
        } else logexit[i] = std::log(exit[i]);
    }
    arma::vec dchaz2entry = rate % logentry % dchaz1entry;
    arma::vec dchaz2exit = rate % logexit % dchaz1exit;

    arma::vec drate = dloghaz1 - dchaz1exit + dchaz1entry;
    arma::vec dshape = dloghaz2 - dchaz2exit + dchaz2entry;
    arma::mat dX = X;
    arma::vec d = drate % rate;
    for (unsigned i = 0; i < X.n_cols; i++) {
        dX.col(i) =  d % X.col(i);
    }
    arma::mat dZ = Z;
    d = dshape % shape;
    for (unsigned i = 0; i < Z.n_cols; i++) {
        dZ.col(i) =  d % Z.col(i);
    }
    arma::mat res = arma::join_horiz(dX, dZ);
    if (indiv) return res;
    return arma::sum(res, 0);
}
