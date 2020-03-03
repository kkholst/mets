#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <complex.h>
#include "twostage.h"

using namespace arma;
using namespace Rcpp;


RcppExport SEXP RVpairs2DIMRV(SEXP irvs, SEXP inrvs) 
{
 IntegerVector nrvs(inrvs);
 mat rvs  = Rcpp::as<mat>(irvs);

 rvs.print(); 

for (unsigned j=0;j<rvs.n_rows;j++) { 
	// 2-dimensional array with (2xrandom effects) written as row 
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	
//	Rprintf(" %d \n",lnrv); 
	mat rv= reshape(rvs.row(j),2,lnrv); 
	rv.print(); 
        Rprintf("==============================\n"); 
}

List res; 
res["mm"]=rvs;

  return(res); 
}

