#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include "twostage.h"

using namespace arma;

RcppExport SEXP Bhat(SEXP ds,
		     SEXP Xs, 
		     SEXP theta, 		     
		     SEXP id, 
		     SEXP ididx,
		     SEXP idsize) { // {{{ 
  try {

    uvec           event = Rcpp::as<uvec>(ds); 
    mat                X = Rcpp::as<mat>(Xs);
    mat               Xc = zeros<mat>(X.n_cols,X.n_cols);
    double      thetahat = Rcpp::as<double>(theta);
    unsigned  stop,start = X.n_rows;
    uvec        eventpos = find(event==1);
    mat               dB = zeros(eventpos.n_elem,X.n_cols);
    uvec         cluster = Rcpp::as<uvec>(id);

    uvec clustersize, clustpos;
    umat clusterindex;
    bool doclust = (Rf_isNull)(idsize);
    if (!doclust)  {
      clustersize        = Rcpp::as<uvec>(idsize);
      clusterindex       = Rcpp::as<umat>(ididx);
    }

    // Obtain usual estimates of increments, dB, of the
    // cumulative time-varying effects in Aalens Additive Model
    for (unsigned ij=0; ij<eventpos.n_elem; ij++) {
      unsigned rij = eventpos.n_rows-ij-1;
      stop  = start-1;
      start = eventpos(rij);
      mat Xij = X(span(start,stop), span::all);
      Xc = Xc + Xij.st()*Xij;
      mat U, V; vec s; 
      svd(U,s,V,Xc);
      mat Xci = U*diagmat(1/s)*V.st();
      dB.row(rij) = trans(Xci*trans(X.row(start)));
    }
    if (thetahat==0) {
      return(Rcpp::List::create(Rcpp::Named("dB")=dB)); //Increments of marg. aalenn
    }


    vec       Hij(X.n_rows); // Vector to hold cumulative hazard; updated as t increases
    Hij.fill(0);
    mat        B2 = zeros(eventpos.n_elem,X.n_cols); // Cond. cumulative coef 
    for (unsigned k=0; k<eventpos.n_elem; k++) { // Iterate over events
      unsigned ij = eventpos(k);
      unsigned i = cluster(ij);  // cluster
      if (doclust) {
         clustpos = find(cluster==i); // position of within cluster observations
      } else {
        unsigned csize = clustersize(i);
        clustpos  = conv_to<uvec>::from(clusterindex.submat(i,0,i,csize-1));
      }
      uvec posL = find(clustpos<ij); // earlier events/censoring within cluster
      uvec posR = find(clustpos>=ij); // later/current events within cluster
      unsigned Ni = 0; // Number of events in cluster before current event,time t-
      double Hi = 0 ; // Sum of cum.haz. within cluster up to time t-
      if (posL.n_elem>0) {
	Ni = sum(event.elem(clustpos.elem(posL)));
	Hi = sum(Hij.elem(clustpos.elem(posL)));
      }
      uvec pos;
      if (posR.n_elem>0 && k>0) {
	pos = clustpos.elem(posR);
	mat Xi = X.rows(pos);
	Hij.elem(pos) = Xi*trans(B2.row(k-1));
	Hi += sum(Hij.elem(pos));
      }
      double psi = (1/thetahat+Ni)/(1/thetahat+fmax(0,Hi));
      B2.row(k) = dB.row(k)/psi;
      if (k>0) { B2.row(k) += B2.row(k-1); }
    }
    return(Rcpp::List::create(Rcpp::Named("dB")=dB, //Increments of marg. aalen
			      Rcpp::Named("B2")=B2  // Cum.coef of frailty aalen
			      ));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

RcppExport SEXP Uhat(SEXP ds, SEXP H, SEXP theta, SEXP id, SEXP idsize) {
  try { // {{{ 
    uvec                event = Rcpp::as<uvec>(ds);
    vec                   Hij = Rcpp::as<vec>(H);
    double           thetahat = Rcpp::as<double>(theta);
    umat              cluster = Rcpp::as<umat>(id);
    uvec clustersize, ucluster, clustpos;
    unsigned nclust;
    bool doclust = (Rf_isNull)(idsize);
    //bool doclust = (cluster.n_cols==1);
    if (doclust)  {
      ucluster              = unique(cluster);
      nclust                = ucluster.n_elem;
    } else {
      clustersize           = Rcpp::as<uvec>(idsize);
      nclust = cluster.n_rows;
    }
    vec res(nclust);
    for (unsigned i=0; i<nclust; i++) {
      if (doclust) {
        unsigned ic = ucluster(i);  // cluster 
        clustpos = find(cluster==ic); // position of within cluster observations
      } else {
        unsigned csize = clustersize(i);
        clustpos  = conv_to<uvec>::from(cluster.submat(i,0,i,csize-1));
	//cluster(span(i,i),span(0,csize-1));
      }
      double Ni = sum(event.elem(clustpos));
      double Hi = sum(Hij.elem(clustpos));
      double thetaH = thetahat*Hi+1;
      double R = (log(thetaH)/thetahat + (Ni-Hi)/(thetaH));
      for (unsigned h=0; h<Ni; h++) R -= 1/(1+thetahat*h);
      res(i) = R/thetahat;
    }
    return(Rcpp::wrap(res));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}} 


RcppExport SEXP MatxCube(
		SEXP imat,
		SEXP idim,SEXP iDBhat 
		)  
{ // {{{ 
try {

 mat  xmat = Rcpp::as<mat>(imat);

 NumericVector vDBhat(iDBhat);
 IntegerVector arrayDims(idim);
 arma::cube DBhat(vDBhat.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

 mat X(arrayDims[2],arrayDims[0]); 

 for (int k=0; k<arrayDims[2]; k++) { // Iterate over events
//	    X.row(k)=  DBhat.slice(k) * trans(xmat.row(k)); 
	    X.row(k)=  xmat.row(k) * trans(DBhat.slice(k)) ;
 } 

    return(Rcpp::List::create(Rcpp::Named("X")=X));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

RcppExport SEXP BhatAddGam(SEXP irecursive, SEXP idBaalen,SEXP icause,
		SEXP idimxjump,SEXP ixjump, // cube 
		SEXP itheta,
		SEXP idimthetades,SEXP ithetades, // cube 
		SEXP iags, SEXP ivarlink, 
		SEXP idimjumprv,SEXP ijumprv,SEXP iit, SEXP iBit)  // cube 
{ // {{{ 
  try {

//   wall_clock timer; 
//   timer.tic(); 

// {{{ reading in matrices and cubes 
    mat                dBaalen = Rcpp::as<mat>(idBaalen);
    uvec                 cause = Rcpp::as<uvec>(icause); 
    vec                  theta = Rcpp::as<vec>(itheta); 
    mat                    ags = Rcpp::as<mat>(iags);
    int                varlink = Rcpp::as<int>(ivarlink);
    int                     it = Rcpp::as<int>(iit);
    int              recursive = Rcpp::as<int>(irecursive);
    mat                    Bit = Rcpp::as<mat>(iBit);

 if (recursive==1) it=1; 

// array for xjump covariates of jump subject, for all causes 
 NumericVector vxjump(ixjump);
 IntegerVector arrayDims(idimxjump);
 arma::cube xjump(vxjump.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vecthetades(ithetades);
 IntegerVector arrayDims1(idimthetades);
 arma::cube thetades(vecthetades.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vrv(ijumprv);
 IntegerVector arrayDims2(idimjumprv);
 arma::cube rv(vrv.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);

 // }}}
 
// double nt = timer.toc();
// printf("timer-ind %lf \n",nt); 

  vec casev(cause.n_elem); 
  vec etheta=theta; 
  if (varlink==1) etheta=exp(theta); 

//  xjump for each jump contains matrix of covariates such that vec cumhaz1= xjump.slice(s) * Bhat 

    mat  Bhat(dBaalen.n_rows, dBaalen.n_cols); 
//    mat  Bhatmarg(dBaalen.n_rows, dBaalen.n_cols); 
//    cube  DthetaBhat(theta.n_elem, dBaalen.n_cols,dBaalen.n_rows); 
//    vec dBB(theta.n_elem); 

//    Bhat.fill(0); // initialize 
//
    vec DthetaS(theta.n_elem),DthetaDtS(theta.n_elem),DthetaW(theta.n_elem); 
    vec allvec(6); 
    int ncr=rv.n_rows; 
    vec cumhaz(ncr); cumhaz.fill(0); 
    vec Dcumhaz1(ncr); 
//    vec cumhaz2(ncr); cumhaz2.fill(0); 
    double  caseweight=1,ll; 
//    mat rv2=0*rv.slice(0); 
      mat rv1=rv.slice(0); 

//	rv1.print("rv1"); 
//	cumhaz1.print("ch1"); 
//	ags.print("ags"); 

// wall_clock timer; 
// timer.tic(); 

    for (int i=0; i<it ; i++) { // Iterate over baseline 
    for (unsigned k=0; k<cause.n_elem; k++) { // Iterate over events
   // {{{ 

        if (recursive==0) cumhaz=xjump.slice(k) * trans(Bit.row(k)); 
        // computes weights based on additive gamma model 
        mat thetadesv=thetades.slice(k); 
	rv1=rv.slice(k); 
        ll=survivalRVCmarg(etheta,thetadesv,ags,(int) cause(k),cumhaz,rv1,DthetaS,DthetaDtS,allvec);
        caseweight=allvec(0)/ll; //   S / D_1 S
	casev(k)=caseweight; 
	vec dbb=trans(dBaalen.row(k)); 
//	Rprintf("%d %d %lf %lf \n",k,cause(k),ll,allvec(0),dbb(0),dbb(1));  
//	etheta.print("theta"); 
//	thetadesv.print("thetades"); 
//	cumhaz.print("cumhaz"); 
//	rv1.print("rv1"); 

//	DthetaW=(ll*DthetaDtS-allvec(0)*DthetaS)/(ll*ll);

        //  increments 
        Bhat.row(k)=dBaalen.row(k)*caseweight;
//        DthetaBhat.slice(k)= DthetaW * dBaalen.row(k);

	// derivative of baseline wrt theta
//	mat dBthetamat=xjump.slice(k) * trans(DthetaBhat.slice(k)); 
//	dBB = trans(dBthetamat) *Dcumhaz1; 

//  cumulative  for all causes
        if (k>0) { Bhat.row(k) += Bhat.row(k-1); }
//        if (k>0) { DthetaBhat.slice(k)+= DthetaBhat.slice(k-1)+ 
//                      (dBB * dBaalen.row(k)); 
//	}
//	mat xj=xjump.slice(k); 
//	xj.print("xj"); 
//	vec bb=trans(Bhat.row(k)); 
//	bb.print("bb"); 
//	cumulative hazard at time t- for all causes 
	if (recursive==1) cumhaz=xjump.slice(k) * trans(Bhat.row(k)); 
	// update Bit Bit.row(k)=Bhat.row(k); 
		
//	mat pp=DthetaBhat.slice(k); 
//	pp.print("pp"); Dcumhaz1.print("Dcumhaz1");  
    } // }}}
    } 

// double nt2 = timer.toc();
// printf("Bhat-profile timer-loop %lf \n",nt2); 

    return(Rcpp::List::create(Rcpp::Named("B")=Bhat, 
			      Rcpp::Named("caseweights")=casev)
//			      Rcpp::Named("DthetaBhat")=-1*DthetaBhat)
		    );
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}


// marginal hazard estimation via iterative estimator
// pairs where we condition on second subjects 
// ascertainment correction is equivalent to case-control sampling 
// (except for delayed entry)
RcppExport SEXP BhatAddGamCC(SEXP itwostage,SEXP idBaalen,SEXP icause,
		SEXP idimxjump,SEXP ixjump, // cube 
		SEXP itheta,
		SEXP idimthetades,SEXP ithetades, // cube 
		SEXP iags, SEXP ivarlink, 
		SEXP idimjumprv,SEXP ijumprv,
		SEXP inrvs, SEXP iit, SEXP iBit, 
		SEXP iBcaseit, SEXP icausecase)  // cube 
{ // {{{ 
  try {

//   wall_clock timer; 
//   timer.tic(); 

// {{{ reading in matrices and cubes 
    int               twostage = Rcpp::as<int>(itwostage);
    mat                dBaalen = Rcpp::as<mat>(idBaalen);
    uvec                 cause = Rcpp::as<uvec>(icause); 
    vec                  theta = Rcpp::as<vec>(itheta); 
    mat                    ags = Rcpp::as<mat>(iags);
    //int                varlink = Rcpp::as<int>(ivarlink);
    uvec                  nrvs = Rcpp::as<uvec>(inrvs); 
    int                     it = Rcpp::as<int>(iit);
    mat                    Bit = Rcpp::as<mat>(iBit);
    mat                Bcaseit = Rcpp::as<mat>(iBcaseit);
    uvec             causecase = Rcpp::as<uvec>(icausecase); 

    // if it>=1 then uses iterative estimator rather than recursive
    // for baseline, case/proband is handled via iterative
    int recursive=1;int ascertainment=1; 
    if (twostage==-2)  { recursive=1; ascertainment=1;it=1;} // not two stage model 
    if (twostage==-1)  { recursive=0; ascertainment=0;     } // not two stage model 
    if (twostage==0)   { recursive=1; ascertainment=0;it=1;} // not two stage model 
    if (twostage==1)   { recursive=1; ascertainment=0;it=1;} // two-stage 
    if (twostage==2)   { recursive=0; ascertainment=0;     } // two-stage 
    if (twostage==3)   { recursive=1; ascertainment=1;it=1;} // two-stage it=1; 

//  xjump for each jump contains matrix of covariates such that vec cumhaz1= xjump.slice(s) * Bhat 
// array for xjump covariates of jump subject, for all causes 
 NumericVector vxjump(ixjump);
 IntegerVector arrayDims(idimxjump);
 arma::cube xjump(vxjump.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vecthetades(ithetades);
 IntegerVector arrayDims1(idimthetades);
 arma::cube thetades(vecthetades.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vrv(ijumprv);
 IntegerVector arrayDims2(idimjumprv);
 arma::cube rv(vrv.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);

 // }}}
 
// double nt = timer.toc(); printf("timer-ind %lf \n",nt); 

  vec casev(cause.n_elem); 
  vec etheta=theta; 
//  if (varlink==1) etheta=exp(theta); 

  mat  Bhat(dBaalen.n_rows, dBaalen.n_cols); 
//    cube  DthetaBhat(theta.n_elem, dBaalen.n_cols,dBaalen.n_rows); 
//    vec dBB(theta.n_elem); 
//    Bhat.fill(0); // initialize 

  vec DthetaS(theta.n_elem); // ,DthetaDtS(theta.n_elem),DthetaW(theta.n_elem); 
  vec allvec(6),  allvect(6); 
  int ncr=rv.n_rows/2; 
//    printf(" %d \n",ncr); 
  vec cumhaz(ncr); cumhaz.fill(0); 
  vec cumhazt(ncr); cumhazt.fill(0); 
//    vec Dcumhaz1(ncr); 
//    vec cumhaz2(ncr); cumhaz2.fill(0); 

  double s1=1,s2=1,caseweight=1,ll=1; 
  double llt=1; //s1t=1,
//    mat rv2=0*rv.slice(0); 
//    mat rv1=rv.slice(0); 

//	rv1.print("rv1"); cumhaz1.print("ch1"); ags.print("ags"); 

// wall_clock timer; 
// timer.tic(); 
 mat rrv1,rrv2; 
 vec rv1,rv2; 

    for (int i=0; i<it ; i++) { // Iterate over baseline  or not
    for (unsigned k=0; k<cause.n_elem; k++) { // Iterate over events
    // {{{ 
    
        // computes weights based on additive gamma model 
	// cumulative hazards for cases Xcase^T B(T_case)
        vec cumhazcase=trans(Bcaseit.row(k)); 
        if ((recursive==0) && (k>0)) cumhaz=xjump.slice(k)*trans(Bit.row(k));
        cumhazt=xjump.slice(k)*trans(Bit.row(k)); 

        //int lnrv= nrvs(k)-1; // number of random effects for this cluster 	
	mat rvm=rv.slice(k);
	int nnn=rvm.n_rows; 
//	rvm.print("rvm"); 
	// first half of rows for person1 and second half for subject 2
	// nn number of competing risks
        if (twostage<=0) {	
           rrv1= rvm.rows(0,nnn/2-1);
           rrv2= rvm.rows(nnn/2,nnn-1);
	} else {
           rv1= trans(rvm.row(0));
           rv2= trans(rvm.row(1));
//	   rv1.print("rv1");  rv2.print("rv2"); 
//	rv1.print("rv1"); thetadesv.print("thetades"); 
	}

        mat thetadesv=thetades.slice(k); 

        if (twostage<=0) {
           ll=survivalRVC2all(etheta,thetadesv,ags,cause(k),causecase(k),cumhaz,cumhazcase,rrv1,rrv2,DthetaS,allvec);
	   if (ascertainment==0) caseweight=allvec(4)/(ll); 
           if (ascertainment==1) {
//	       s1t=exp(-cumhazt(0)); 
	       llt=survivalRVC2all(etheta,thetadesv,ags,0,causecase(k),cumhaz,cumhazcase,rrv1,rrv2,DthetaS,allvect);
	       caseweight=(allvec(4)+llt)/(ll); 
	    }
	} else { 
	    s1=exp(-cumhaz(0));      // survival to clayton-oakes
	    s2=exp(-cumhazcase(0));  // cases via iterative, one-dimensional only (survival)
//         printf(" %d %d %lf %lf \n",cause(k),causecase(k),s1,s2); 
//	    etheta.print("hlj");  ags.print("ags"); 
	    ll=claytonoakesRVC(etheta,thetadesv,ags,cause(k),causecase(k),s1,s2,rv1,rv2,DthetaS,allvec);
	    if (ascertainment==0) caseweight=allvec(0)/(s1*ll); 
	    if (ascertainment==1) {
//	       s1t=exp(-cumhazt(0)); 
	       llt=claytonoakesRVC(etheta,thetadesv,ags,0,causecase(k),s1,s2,rv1,rv2,DthetaS,allvect);
	       caseweight=(s1*allvec(0)+s2*llt)/(s2*s1*ll); 
	    }
//	    printf("caseweight %lf %lf \n",caseweight,ll); 
	}
	//   either S / D1 S, when cause2=0, or D2 S / D1D2S, when cause2=1
	casev(k)=caseweight; 

        //  increments 
        Bhat.row(k)=dBaalen.row(k)*caseweight;
//        DthetaBhat.slice(k)= DthetaW * dBaalen.row(k);

	// derivative of baseline wrt theta
//	mat dBthetamat=xjump.slice(k) * trans(DthetaBhat.slice(k)); 
//	dBB = trans(dBthetamat) *Dcumhaz1; 

//  cumulative  for all causes
        if (k>0) { Bhat.row(k) += Bhat.row(k-1); }
//      if (k>0) { DthetaBhat.slice(k)+= DthetaBhat.slice(k-1)+ 
//                      (dBB * dBaalen.row(k)); 
//	}
	mat xj=xjump.slice(k); 
//	vec bb=trans(Bhat.row(k)); 
//	bb.print("bb"); 
//	cumulative hazard at time t- for all causes 
	if (recursive==1)  cumhaz=xjump.slice(k)*trans(Bhat.row(k)); 
//	else  { 
//		printf("sono qui \n"); 
		// update Bit Bit.row(k)=Bhat.row(k); 
//	}
		
//	mat pp=DthetaBhat.slice(k); 
//	pp.print("pp"); Dcumhaz1.print("Dcumhaz1");  

//		trans(sum(dBthetamat,0)); 
    } // }}}
    } 

// double nt2 = timer.toc();
// printf("Bhat-profile timer-loop %lf \n",nt2); 

    return(Rcpp::List::create(Rcpp::Named("B")=Bhat, 
			      Rcpp::Named("caseweights")=casev)
//			      Rcpp::Named("DthetaBhat")=-1*DthetaBhat)
		    );
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

RcppExport SEXP XBmindex(SEXP imindex,SEXP iX,SEXP icumt)  
{ // {{{ 
try {
 mat  index = Rcpp::as<mat>(imindex);
 mat  cumt  = Rcpp::as<mat>(icumt);
 mat  X     = Rcpp::as<mat>(iX);
 int n=index.n_rows; 
 int p=X.n_cols; 

 mat XBmindex(n,n); 
 vec vcumt(p); 

// XBmindex.print("XB"); 

 for (int r=0; r<n; r++) { // Iterate over events
   rowvec Xvec=X.row(r); 
 for (int c=0; c<n; c++) { // Iterate over events
   int ii=index(r,c)-1; 
   if (ii>0) {
   vcumt=trans(cumt.row(ii)); 
//   vcumt.print("vcum"); 
//   Xvec.print("Xv"); 
//   printf(" %d %d \n",r,c); 
   mat  out = Xvec*vcumt; 
   XBmindex(r,c)= out(0,0); 
   }
 } 
 }

    return(Rcpp::List::create(Rcpp::Named("XBmindex")=XBmindex));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

//RcppExport SEXP backfitEaEt(SEXP iYt,SEXP iage0, SEXP iagenull, SEXP iage, 
//		SEXP itime0, SEXP itime, SEXP iajumps, SEXP itjumps,SEXP ia0)  
//{ // {{{ 
//try {
//
// vec  Yt  = Rcpp::as<vec>(iYt);
// vec  age = Rcpp::as<vec>(iage); 
// vec  age0 = Rcpp::as<vec>(iage0); 
// vec  agenull = Rcpp::as<vec>(iagenull); 
// vec  time = Rcpp::as<vec>(itime); 
// vec  time0 = Rcpp::as<vec>(itime0); 
// vec  ajumps = Rcpp::as<vec>(iajumps); 
// vec  tjumps = Rcpp::as<vec>(itjumps); 
// double a0= Rcpp::as<double>(ia0);
// 
// int n=Yt.n_rows; 
// mat Ea(n,n); Ea.zeros(); 
// mat Et(n,n); Et.zeros(); 
// double ao,to,nn,ttt; 
//
// for (int i=0; i<n; i++) { 
// for (int j=0; j<n; j++) { // Iterate over events
//     ao=ajumps(i); to=tjumps(j);
//     vec who=(to+agenull>age0)*(to+agenull<=age); 
//     vec top=who*(a0-agenull<to)*(to<=ao-agenull); 
//     nn=sum(who); 
//     ttt=sum(top)/nn; 
//     if (nn>0) Et(i,j)=ttt; 
//     vec aaa=ao-agenull; 
//// ###     aaa[aaa<0] <- 0
////     arma::uvec fff = FastApproxC(jumps,aaa,TRUE,2); 
////     vec yaai=Yt.elem(fff); 
//     vec yaai=(top>0);  
//     vec ww=(yaai>0); 
////     ### at risk at time a
//     who=(ao-agenull>time0)*(ao-agenull<time); 
//     top= who*(agenull<ao)*(ao<=agenull+to); 
//     ttt=0; 
//     for (int k=0;k<n; k++) 
//	     if (yaai(k)>0) ttt+=top(k)/yaai(k); 
//     Ea(j,i) = ttt; 
// }
// }
//
//    return(Rcpp::List::create(Rcpp::Named("Et")=Et, Rcpp::Named("Ea")=Ea));
//  } catch( std::exception &ex ) {
//    forward_exception_to_r( ex );
//  } catch(...) {  
//    ::Rf_error( "c++ exception (unknown reason)" ); 
//  }
//  return R_NilValue; // -Wall
//} // }}}
