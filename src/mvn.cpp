// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "mvn.h"
#include "tools.h"
#include <math.h>

#include <R_ext/Rdynload.h>  /* required by R */
#include <mvtnormAPI.h>
#include <mvtnorm.h>

int _mvt_maxpts=25000;
double _mvt_abseps=0.001;
double _mvt_releps=0;
int _mvt_df = 0;
int _mvt_inform;
double _mvt_error[3];


double mvtdst(int* n,
	      int* nu,        // Degrees of freedom (0=MVN)
	      double* lower,  // Lower integration bounds
	      double* upper,  // Upper integration bounds
	      int* infin,     // Infinity argument, ith element 0: ]-inf,up[i]], 1: [lo[i],inf[, 2: [lo[i],up[i]], 3: ]-inf,inf[
	      double* correl, // Correlation coefficients (upper-tri)
	      double* delta,  // non-central parameter
	      int* maxpts,    // Max function evalutions (quasi-mc)
	      double* abseps, // Tolerance absolute error
	      double* releps, // Tolerance relative error
	      double* error,  // estimated abs. error. with 99% confidence interval
	      double* value,  // result
	      int* inform)    // Message (0 success)
{
  if (*n==1 && *nu==0) {
    // 0: right, 1: left, 2: interval
    switch (*infin) {
    case 0: *value = Rf_pnorm5(*upper,0.0,1.0,1,0); break;
    case 1: *value = 1-Rf_pnorm5(*lower,0.0,1.0,1,0); break;
    case 2: *value = Rf_pnorm5(*upper,0.0,1.0,1,0)-Rf_pnorm5(*lower,0.0,1,1,0); break;
      //    default: *value = 0;
    }
    //cerr << "infin=" << *infin << " value=" << *value << " upper=" << *upper << " lower=" << *lower << endl;
    return(*value);
  }

    int rnd=1;
  /* mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h */
  mvtnorm_C_mvtdst(n, nu,
		   lower, upper, infin, correl, delta,
		   maxpts, abseps, releps,
		   error, value, inform, &rnd);
  switch (*inform) {
  case 0:
    return *value;
  case 1:
  case 2:
  case 3:
    return -1.0;
  };
  return *value;
}


double dmvn(const vec &y, const mat &W,
	    bool log=true, double logdet=datum::inf) {
  int n = W.n_rows;
  double res = -0.5*n*log2pi;
  if (logdet!=datum::inf) {
    res += -0.5*(logdet + as_scalar(trans(y)*W*y));
  } else {
    double sign=0;
    mat iW = inv(W);
    log_det(logdet,sign,W);
    res += -0.5*(logdet + as_scalar(trans(y)*iW*y));
  }
  if (!log) res = exp(res);
  return(res);
}


double cdfmvn(mat &upper, mat &cor) {
  double val=0;
  int n = cor.n_cols;
  rowvec _mvt_delta(n); _mvt_delta.fill(0);
  unsigned ncor = n*(n-1)/2;
  rowvec Cor(ncor);
  int j = 0;
  for (int r=0; r<n; r++) {
    for (int c=r+1; c<n; c++) {
      Cor(j) = cor(r,c);
      j++;
    }
  }
  Row<int> infin(n); infin.fill(0); //
  mvtdst(&n,
	 &_mvt_df,
	 &upper[0], // Lower, ignored
	 &upper[0],
	 &infin[0], // Infinity argument (all 0 since CDF)
	 &Cor[0],
	 &_mvt_delta[0],
	 &_mvt_maxpts,
	 &_mvt_abseps, &_mvt_releps,
	 &_mvt_error[0],
	 &val, &_mvt_inform);
  return(val);
}



RcppExport SEXP Dpmvn (SEXP lower, SEXP upper, SEXP mu, SEXP sigma, SEXP std) {
BEGIN_RCPP
  colvec y = Rcpp::as<colvec>(upper);
 NumericVector Lower(lower);
 bool Std = Rcpp::as<bool>(std);
 mat S = Rcpp::as<mat>(sigma);
 colvec Mu = Rcpp::as<colvec>(mu);
 int n = S.n_cols;
 vec L(n);
 for (int j=0; j<n; j++) {
   if (y(j)==datum::inf) L(j) = 1; // (lower,+Inf)
   if (Lower(j)==-datum::inf) L(j) = 0; // (-Inf,upper)
   // =2 (low,up)   For now we actually don't consider the interval-censored case
    // <0 (-Inf,Inf)
 }
 mat LL = diagmat(L);
 mat iLambda;
 if (!Std) {
   iLambda = diagmat(1/sqrt(diagvec(S)));
   S = iLambda*S*iLambda;
   y = iLambda*(y-Mu);
 }

 List res;

 vec D(n);
 mat H(n,n);
  if (n==1) {
    List res;
    D[0] = Rf_dnorm4(y[0],0.0,1.0,0);
    H[0] = -y[0]*D[0];
    if (!Std) {
      D[0] *= iLambda[0];
      H[0] *= iLambda[0]*iLambda[0];
    }
    res["gradient"] = D;
    res["hessian"] = H;
    return(res);
  }

  for (int j=0; j<n; j++) {
    mat Sj = S; Sj.shed_row(j);
    mat Sj0 = Sj.col(j);
    Sj.shed_col(j);
    Sj -= Sj0*Sj0.t();
    mat iL = diagmat(1/sqrt(diagvec(Sj)));
    Sj = iL*Sj*iL;
    vec muj = y;
    muj.shed_row(j);
    muj = iL*(muj - Sj0*y[j]);
    D[j] = Rf_dnorm4(y[j],0.0,1.0,0)*
      cdfmvn(muj,Sj);
  }

  if (n==2) {
    H(0,1) = H(1,0) = dmvn(y, S, false);
    H.diag() = -y%D - H(0,1)*S(0,1);
  } else {
    mat phis(n,n); phis.fill(0);
    mat Phis = phis;
    uvec idx1(n-2);
    uvec idx2(2);
    for (int i=0; i<(n-1); i++) {
      for (int j=(i+1); j<n; j++) {
	idx2(0) = i; idx2(1) = j;
	unsigned pos = 0;
	for (int k=0; k<n; k++) {
	  if (k!=i && k!=j) {
	    idx1(pos) = k;
	    pos++;
	  }
	}
	mat Snij = S.submat(idx1,idx2);
	mat S2 = S.submat(idx2,idx2);
	mat B = Snij*inv(S2);
	mat Sij = S.submat(idx1,idx1) - B*trans(Snij);
	vec muij = y.elem(idx1) - B*y.elem(idx2);

	mat iL = diagmat(1/sqrt(diagvec(Sij)));
	Sij = iL*Sij*iL;
	muij = iL*muij;
	Phis(i,j) = Phis(j,i) = cdfmvn(muij,Sij);
	phis(i,j) = phis(j,i) = dmvn(y.elem(idx2),S2,false);
      }
    }
    H = Phis%phis;
    colvec ones(n); ones.fill(1);
    H.diag() = -y%D - (S%H)*ones;
  }

  //double val = cdfmvn(y,S);
  //return(Rcpp::wrap(val));
  if (!Std) {
    H = iLambda*H*iLambda;
    D = iLambda*D;
  }
  res["gradient"] = D;
  res["hessian"] = H;

  return(res);
END_RCPP
}



// RcppExport SEXP scoremvn2 (SEXP lower, SEXP upper, SEXP mu, SEXP sigma, SEXP std) {
// BEGIN_RCPP
  // int k = Yl.n_cols;
  // int n = Yl.n_rows;
  // uvec Obs = find(Status==0);
  // uvec Cens = find(Status==1);
  // uvec Ord = find(Status>1);
  // uvec NonObs = find(Status>0);
  // int nObs = Obs.size();
  // int nNonObs = NonObs.size();
  // int nOrd = Ord.size();
  // int nu = Su.n_cols;
  // bool nonconstvar = (nu>0);

  // double sign, logdetS0;
  // mat Se,S0,iS0;

  // vec loglik(n); loglik.fill(0);
  // //  return(loglik);
  // if (nObs>0) {
  //   mat Y0 = Yl.cols(Obs)-Mu.cols(Obs);
  //   Se = S0 = S.submat(Obs,Obs);
  //   iS0 = inv(S0);
  //   log_det(logdetS0,sign,S0);
  //   double normconst = -0.5*nObs*log2pi;

  //   for (int i=0; i<n; i++) { // Iterate over subjects

  //     if (nonconstvar) {
  // 	mat Z0 = reshape(Z.row(i),k,nu);
  // 	S0 = Se+Z0.rows(Obs)*Su*trans(Z0.rows(Obs));
  // 	iS0 = inv(S0);
  // 	log_det(logdetS0,sign,S0);
  //     }

  //     loglik(i) = -0.5*(logdetS0 + as_scalar(Y0.row(i)*iS0*trans(Y0.row(i))));
  //     loglik(i) += normconst;
  //   }
//   }

// END_RCPP


// [[Rcpp::export(name = ".scoreMVN")]]
arma::mat scoreMVN(arma::mat &Y,
                   arma::mat &Mu,
                   arma::mat &dMu,
                   arma::mat &S,
                   arma::mat &dS,
                   double itol=0.0) {
  // mat &Z,  mat &Su, mat &dSu,
  // mat &Threshold, mat &dThreshold) {

  int n = Y.n_rows;
  int p = dMu.n_cols;
  double logdet = 0;

  mat iS = Inv(S, logdet, itol);
  mat z = (Y-Mu)*iS;
  mat U(n,p);
  U.fill(0);
  colvec A = -0.5*trans(dS)*vectorise(iS);
  for (int i=0; i<n; i++) {
    colvec zi = trans(z.row(i));
    U.row(i) = trans( A+0.5*trans(dS)*vectorise(zi*trans(zi))
                      + trans(dMu)*zi ) ;
  }
  return(U);
}

vec loglikmvn(mat &Yl, mat &Yu, uvec &Status, mat &Mu, mat &S,
              mat &Z, mat &Su,
              mat &Threshold,
              double itol=0.0, bool nonconstvar = false) {

  int k = Yl.n_cols;
  int n = Yl.n_rows;
  uvec Obs = find(Status==0);
  uvec Cens = find(Status==1);
  uvec Ord = find(Status>1);
  uvec NonObs = find(Status>0);
  int nObs = Obs.size();
  int nNonObs = NonObs.size();
  int nOrd = Ord.size();
  double logdetS0;
  mat Se,S0,iS0;

  vec loglik(n); loglik.fill(0);
  //  return(loglik);
  if (nObs>0) {
    mat Y0 = Yl.cols(Obs)-Mu.cols(Obs);
    Se = S0 = S.submat(Obs,Obs);
    iS0 = Inv(S0, logdetS0, itol);
    // iS0 = inv(S0);
    // log_det(logdetS0, sign, S0);
    double normconst = -0.5*nObs*log2pi;

    for (int i=0; i<n; i++) { // Iterate over subjects

      if (nonconstvar) {
        unsigned nu = Su.n_cols;
        mat Z0 = reshape(Z.row(i),k,nu);
        S0 = Se+Z0.rows(Obs)*Su*trans(Z0.rows(Obs));
        iS0 = Inv(S0, logdetS0, itol);
      }

      loglik(i) = -0.5*(logdetS0 + as_scalar(Y0.row(i)*iS0*trans(Y0.row(i))));
      loglik(i) += normconst;
    }
  }


  if (nNonObs>0) {
    rowvec _mvt_delta(nNonObs); _mvt_delta.fill(0);

    mat MuNonObs = Mu.cols(NonObs);
    mat SNonObs = S.submat(NonObs,NonObs);

    if ((nObs>0) & (!nonconstvar)) { // Calculate conditional on observed
      mat S01 = S.submat(NonObs,Obs);
      mat iS1 = Inv(S.submat(Obs,Obs),logdetS0, itol);
      // mat iS1 = inv(S.submat(Obs,Obs));
      MuNonObs = MuNonObs +
        trans(S01*iS1*trans(Yl.cols(Obs)-Mu.cols(Obs)));
      //      MuNonObs.each_row() += Mu.(NonObs);
      SNonObs = SNonObs - S01*iS1*trans(S01);

    }
    vec il = 1/sqrt(diagvec(SNonObs));
    mat iL = diagmat(il);

    Se = S0 = iL*SNonObs*iL; // Correlation matrix
    int ncor = nNonObs*(nNonObs-1)/2;
    rowvec Cor(std::min(1, ncor)); // We allocate a 1x1 matrix even when ncor=0
    if (ncor>0) {
      int j = 0;
      for (int r=0; r<nNonObs; r++) {
        for (int c=r+1; c<nNonObs; c++) {
          Cor(j) = S0(r,c);
          j++;
        }
      }
    }
    int nthresmax = Threshold.n_cols;
    uvec OrdNonObs;

    mat Thres;
    if (nOrd>0) {
      OrdNonObs = find(Status.elem(NonObs)>1);
      Thres = Threshold.rows(Ord);
    }

    rowvec lower(nNonObs);
    rowvec upper(nNonObs);
    Row<int> infin(nNonObs); // 0: right, 1: left, 2: interval

    uvec currow(1);
    for (int i=0; i<n; i++) { // Iterate over subjects

      currow(0) = i;
      rowvec Mi = MuNonObs.row(i);
      lower = Yl.submat(currow,NonObs);
      upper = Yu.submat(currow,NonObs);

      if (nonconstvar) {
        unsigned nu = Su.n_cols;
        mat Z0 = reshape(Z.row(i),k,nu);
        mat SS = S+Z0*Su*trans(Z0);

        if (nObs>0) {
          mat S0 =  SS.submat(NonObs,NonObs);
          mat S01 = SS.submat(NonObs,Obs);
          mat iS1 = Inv(SS.submat(Obs,Obs), logdetS0, itol);
          // mat iS1 = inv(SS.submat(Obs,Obs));
          Mi = Mi +
            trans(S01*iS1*trans(Yl.submat(currow,Obs)-Mu.submat(currow,Obs)));

          SNonObs = S0 - S01*iS1*trans(S01);
        } else {
          SNonObs = SS;
        }

        il = 1/sqrt(diagvec(SNonObs));
        iL = diagmat(il);
        Se = S0 = iL*SNonObs*iL; // Correlation matrix
        if (ncor>0) {
          int j = 0;
          for (int r=0; r<nNonObs; r++) {
            for (int c=r+1; c<nNonObs; c++) {
              Cor(j) = S0(r,c);
              j++;
            }
          }
        }
      }

      umat StatusNonObs = Status.elem(NonObs);
      infin.fill(2);

      for (int j=0; j<nNonObs; j++) {
        if (upper(j)==datum::inf && StatusNonObs(j)==1) infin(j) = 1;
        if (lower(j)==-datum::inf && StatusNonObs(j)==1) infin(j) = 0;
      }
      // uvec infplus = find(upper==datum::inf);
      // uvec infminus = find(lower==-datum::inf);
      // if (infplus.size()>0) { infin.elem(infplus) -= 1; }
      // if (infminus.size()>0) { infin.elem(infminus) -= 2; }

      if (nOrd>0) {
        for (int j=0; j<nOrd; j++) {
          int jj = OrdNonObs(j);
          int yval = (int) lower(jj);
          if (yval<1) { // if Y=0
            infin(jj) = 0; // Integrate over left tail
            upper(jj) = Thres(j,0);
          } else { // Y>1
            if (yval>=nthresmax) { // Y=k (last)
              double val = Thres(j,yval-1);
              infin(jj) = 1; // Integrate over right tail
              lower(jj) = val;
            } else {
              double val = Thres(j,yval-1);
              double val2 = Thres(j,yval);
              if (val>=val2) { // Also Y=k (last)
                infin(jj) = 1;
                lower(jj) = val;
              } else { //Y=k-i (in between)
                lower(jj) = val;
                upper(jj) = val2;
              }
            }
          }
        }
      }

      lower = (lower-Mi)%trans(il);
      upper = (upper-Mi)%trans(il);

      double val;
      val = mvtdst(&nNonObs, &_mvt_df,
                   &lower[0], &upper[0],
                   &infin[0], &Cor[0],
                   &_mvt_delta[0], &_mvt_maxpts,
                   &_mvt_abseps, &_mvt_releps,
                   &_mvt_error[0], &val, &_mvt_inform);

      // if (isnan(val)) {
      // 	cerr << "***i=" << i << endl;
      // 	cerr << "Threshold=" << Threshold << endl;
      // 	cerr << "Thres=" << Thres << endl;
      // 	cerr << "status=" << Status;
      // 	cerr << "yl=" << Yl.row(i);
      // 	cerr << "yu=" << Yu.row(i);
      // 	cerr << "lower=" << lower;
      // 	cerr << "upper=" << upper;
      // 	cerr << "infin=" << infin;
      // 	cerr << "Cor=\n" << Cor;
      // cerr << "   val=" << val << endl;
      // // }
      // cerr << "   1:loglik(i)=" << loglik(i) << endl;
      loglik(i) += log(val);
      // cerr << "   2:loglik(i)=" << loglik(i) << endl;
      // }

    }
  }

  return(loglik);
}


// [[Rcpp::export(name = ".loglikMVN")]]
arma::mat loglikMVN(arma::mat Yl, arma::mat Yu,
                    arma::uvec Status,
                    arma::mat Mu,
                    arma::mat S,
                    arma::mat Threshold,
                    SEXP z, SEXP su,
                    double itol) {

  if ((Mu.n_cols!=Yl.n_cols) || (Mu.n_rows!=Yl.n_rows))
    throw(Rcpp::exception("Dimension of 'mu' and 'yl' did not agree","mvn.cpp",1));
  if (Status.size()!=Yl.n_cols)
    throw(Rcpp::exception("Dimension of 'status' and 'yl' did not agree","mvn.cpp",1));

  uvec Cens = find(Status==1);
  uvec Obs = find(Status==0);
  uvec NonObs = find(Status>0);
  uvec Ord = find(Status>1);
  // int nObs = Obs.size();
  // int nNonObs = NonObs.size();
  int nCens = Cens.size();
  unsigned n = Yl.n_rows;
  mat Z,Zsub;
  mat Su;
  if (!Rf_isNull(z)) {
    Z = Rcpp::as<mat>(z);
    Su = Rcpp::as<mat>(su);
    if (Z.n_cols!=(Yl.n_cols/Su.n_cols))
      throw(Rcpp::exception("Dimension of 'z' and 'su' did not agree","mvn.cpp",1));
    if (Z.n_rows!=n)
      throw(Rcpp::exception("Dimension of 'z' and 'yl' did not agree","mvn.cpp",1));
  }
  // int nOrd = Ord.size();
  // mat Threshold;
  // if (nOrd>0) {
  //   Threshold = Rcpp::as<mat>(threshold);
  //   if (Threshold.n_rows!=Yl.n_cols)
  //     throw(Rcpp::exception("Dimension of 'threshold' and 'yl' did not agree","mvn.cpp",1));
  // }

  vec loglik(n); loglik.fill(0);
  if (nCens>0) {
    if ((Yl.n_cols!=Yu.n_cols) || (Yl.n_rows!=Yu.n_rows))
      throw(Rcpp::exception("Dimension of 'yl' and 'yu' did not agree","mvn.cpp",1));

    umat stat = (Yl.cols(Cens)==Yu.cols(Cens));
    uvec group(n);
    umat pattern;
    fastpattern(stat,pattern,group);
    uvec NewStatus = Status;
    unsigned K = pattern.n_rows;
    for (unsigned i=0; i<K; i++) { // Iterate over patterns
      uvec idx = find(group==i);
      NewStatus.elem(Cens) = 1-pattern.row(i);
      mat Ylsub = Yl.rows(idx);
      mat Yusub = Yu.rows(idx);
      mat Musub = Mu.rows(idx);
      if (!Rf_isNull(z)) Zsub = Z.rows(idx);

      vec ll = loglikmvn(Ylsub, Yusub, NewStatus, Musub, S, Zsub, S, Threshold,
			               itol);
      loglik.elem(idx) = ll;
    }
  } else {
    loglik = loglikmvn(Yl,Yu,Status,
    		       Mu, S,
    		       Z,  Su,
    		       Threshold,
		       itol);
  }

  return (loglik);
}
//   return(Rcpp::List::create(
// 			    Rcpp::Named("loglik")=loglik,
// 			    Rcpp::Named("k")=k,
// 			    Rcpp::Named("n")=n
// 			    ));
// }


void cov2cor0(const mat &x, rowvec &Cor, rowvec &sx, bool nrm=true) {
  unsigned p = x.n_cols;
  if (nrm) for (unsigned j=0; j<p; j++) sx(j) = 1/sqrt(x(j,j));
  unsigned j=0;
  for (unsigned r=0; r<p; r++) {
    for (unsigned c=r+1; c<p; c++) {
      if (nrm)
        Cor(j) = x(r,c)*sx(r)*sx(c);
      else
        Cor(j) = x(r,c);
      j++;
    }
  }
}

RcppExport SEXP pmvn0(SEXP lower, SEXP upper,
		      SEXP mu, SEXP sigma, SEXP cor) {
BEGIN_RCPP
  mat Sigma = Rcpp::as<mat>(sigma);
  bool asCor = Rcpp::as<bool>(cor);
  mat Mu = Rcpp::as<mat>(mu);
  mat Lower = Rcpp::as<mat>(lower);
  mat Upper = Rcpp::as<mat>(upper);
  unsigned n = Mu.n_rows;
  int p = Mu.n_cols;
  unsigned ncor = p*(p-1)/2;
  bool nSigma = false;

  if ((asCor && Sigma.n_rows>1) || (!asCor && Sigma.n_cols==unsigned(p*p) && Sigma.n_rows>1)) {
    nSigma = true;
    n = Sigma.n_rows;
  }

  rowvec _mvt_delta(p); _mvt_delta.fill(0); // Non-centrality par.
  rowvec Cor(ncor); // Vector of correlation coefficients (upper-tri, colwise)
  rowvec L(p); // Std.deviations

  //bool nSigma = Sigma.n_rows==n && Sigma.n_cols!=(unsigned)p;
  if (!nSigma) {
    if (!asCor) {
      cov2cor0(Sigma,Cor,L,true);
    } else {
      Cor = Sigma.row(0);
    }
  }

  Row<int> infin(p); infin.fill(2); //
  for (unsigned j=0; j<(unsigned)p; j++) {
    if (Upper(0,j)==datum::inf) infin(j) = 1;
    if (Lower(0,j)==-datum::inf) infin(j) = 0;
  }

  rowvec Lower0(p);
  rowvec Upper0(p);
  rowvec Mu0 = Mu.row(0);
  vec res(n);
  for  (unsigned i=0; i<n; i++) {
    double val;
    if (Mu.n_rows>1) Mu0 = Mu.row(i);
    if (Lower.n_rows==n) {
      Lower0 = Lower.row(i)-Mu0;
      Upper0 = Upper.row(i)-Mu0;
      infin.fill(2); // (a,b)
      for (unsigned j=0; j<(unsigned)p; j++) {
	if (Upper0(0,j)==datum::inf) infin(j) = 1; // (a,Inf)
	if (Lower0(0,j)==-datum::inf) {
	  if (infin(j)==1) infin(j) = -1; // (-Inf,Inf)
	  else infin(j) = 0; // (-Inf,b)
	}
      }
    } else {
      Lower0 = Lower.row(0)-Mu0;
      Upper0 = Upper.row(0)-Mu0;
    }
    // We use that Phi(a,b,S,mu) = Phi(L(a-mu),L(b-mu),R,0); R=LSL
    if (nSigma) {
      if (asCor) {
	Cor = Sigma.row(i);
      } else { // p*p row
	mat Sigma0 = Sigma.row(i);
	Sigma0.reshape(p,p);
	cov2cor0(Sigma0,Cor,L,true);
      }
    }
    if (!asCor) {
      Lower0 = Lower0%L;
      Upper0 = Upper0%L;
    }
    // std::cerr << "Lower" << Lower0;
    // std::cerr << "Upper" << Upper0;
    // std::cerr << "Infin" << infin;
    // std::cerr << "Cor" << Cor;
    // std::cerr << "mvtdelta" << _mvt_delta;
    // std::cerr << "mvt_df" << _mvt_df;
    mvtdst(&p, &_mvt_df,
	   &Lower0[0], &Upper0[0],
	   &infin[0], &Cor[0],
	   &_mvt_delta[0], &_mvt_maxpts,
	   &_mvt_abseps, &_mvt_releps,
	   &_mvt_error[0], &val, &_mvt_inform);
    res(i) = val;
  }
  return(Rcpp::wrap(res));
END_RCPP
}


//////////////////////////////////////////////////
// Bivariate case
//////////////////////////////////////////////////

RcppExport SEXP bvncdf(SEXP a, SEXP b, SEXP r) {
  double u1 = -Rcpp::as<double>(a);
  double u2 = -Rcpp::as<double>(b);
  double rho = Rcpp::as<double>(r);
  double val = bvnd_(&u1, &u2, &rho);
  NumericVector res(1); res[0] = val;
  return(res);
}

double dbvnorm(double y1, double y2, double R) {
  double detR = 1-R*R;
  // inv(R) = [1 -r; -r 1]/detR (prove by gauss elim.)
  double res = 1/(2*M_PI*sqrt(detR))*exp(-0.5/detR*(y1*y1+y2*y2-2*R*y1*y2));
  return(res);
}

vecmat Dbvn(double y1, double y2, double R) {
  vec DP(2);
  double R2 = R*R;
  DP(0) = Rf_dnorm4(y1,0.0,1.0,0)*Rf_pnorm5(y2-R*y1,0.0,sqrt(1-R2),1,0);
  DP(1) = Rf_dnorm4(y2,0.0,1.0,0)*Rf_pnorm5(y1-R*y2,0.0,sqrt(1-R2),1,0);
  mat HP(2,2);
  HP(1,0) = HP(0,1) = dbvnorm(y1,y2,R);
  HP(0,0) = -y1*DP(0) - R*HP(1,0);
  HP(1,1) = -y2*DP(1) - R*HP(1,0);
  vecmat res;
  res.V = DP;
  res.M= HP;
  return(res);
}

double Sbvn(double &l1, double &l2, double &r) {
  double val = bvnd_(&l1, &l2, &r);
  return(val);
}


//////////////////////////////////////////////////
// Density with varying correlation structure
//////////////////////////////////////////////////

// [[Rcpp::export(name = ".dmvn")]]
NumericVector dmvn(arma::mat u, arma::mat mu, arma::mat rho) {
  unsigned n = u.n_rows;
  unsigned p = u.n_cols;
  NumericVector res(n);
  arma::mat R = arma::eye(p,p);
  double logdetR = 0;
  arma::mat invR = R;
  arma::rowvec mu_i(p);

  for (unsigned i=0; i<n; i++) {
    if (i % 10000 == 0) Rcpp::checkUserInterrupt();
    if (i<mu.n_rows) {
      mu_i = mu.row(i);
    }
    arma::rowvec ui = u.row(i)-mu_i;
    if (i<rho.n_rows) {
      unsigned pos=0;
      for (unsigned r=0; r<p; r++)
	for (unsigned c=r+1; c<p; c++) {
	  R(r,c) = rho(i,pos);
	  R(c,r) = rho(i,pos);
	  pos++;
	}
      logdetR = log(fabs(arma::det(R)));
      invR = inv(R);
    }

    res[i] = -0.5*(logdetR + p*log2pi +
		   arma::as_scalar(ui*(invR)*ui.t()));
  }
  return res;
}

//////////////////////////////////////////////////
// RNG with varying correlation structure
//////////////////////////////////////////////////

// [[Rcpp::export(name = ".rmvn")]]
arma::mat rmvn(unsigned n, arma::mat mu, arma::mat rho) {
  unsigned p = mu.n_cols;
  arma::mat res(n,p);
  std::generate( res.begin(), res.end(), ::norm_rand ) ;
  arma::mat R = arma::eye(p,p);
  arma::mat L = R;
  arma::rowvec mu_i(p);
  arma::rowvec U(p);
  for (unsigned i=0; i<n; i++) {
    if (i % 10000 == 0) Rcpp::checkUserInterrupt();
    if (i<rho.n_rows) {
      unsigned pos=0;
      for (unsigned r=0; r<p; r++)
	for (unsigned c=r+1; c<p; c++) {
	  R(r,c) = rho(i,pos);
	  R(c,r) = rho(i,pos);
	  pos++;
	}
      L = chol(R);
    }
    if (i<mu.n_rows) {
      mu_i = mu.row(i);
    }
    res.row(i) =  res.row(i)*L;
    res.row(i) += mu_i;
  }
  return res;
}
