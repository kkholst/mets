#ifndef TOOLS_H
#define TOOLS_H

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstring>

using namespace Rcpp;

RcppExport SEXP FastLong2(SEXP idata, SEXP inclust, SEXP infixed, SEXP invarying);
RcppExport SEXP FastLong(SEXP idata, SEXP inclust, SEXP infixed, SEXP invarying, SEXP missing);
RcppExport SEXP FastApprox(const SEXP time, const SEXP newtime, const SEXP equal, const SEXP type);
RcppExport SEXP FastPattern(SEXP y1,SEXP y2, SEXP cat);
RcppExport SEXP FastCluster(SEXP x);

void fastpattern(const arma::umat &y, arma::umat &pattern, arma::uvec &group, unsigned categories=2);

template <class T>
std::string numStr(T x) {
  std::ostringstream nmbstr; nmbstr << x;
  std::string ss = nmbstr.str();
  return ss;
}

arma::mat Inv(const arma::mat &X, double &logdet, double itol=0.0);


#endif /* TOOLS_H */
