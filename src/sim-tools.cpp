// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo/Lighter>
#include <Rmath.h>
#include <vector>
using namespace arma; 
using namespace Rcpp; 
// [[Rcpp::depends(RcppArmadillo)]]

 arma::uvec approx(const arma::mat &time,  // Sorted time points (Ascending)
		  const arma::mat &newtime,
		  unsigned type) { // (0: nearest, 1: right, 2: left){{{
    uvec idx(newtime.n_elem);
    double vmax = time(time.n_elem-1);
    vec::const_iterator it;
    double upper=0.0; int pos=0;
    for (int i=0; i<newtime.n_elem; i++) {
      if (newtime[i]>=vmax) {
      pos = time.n_elem-1;
      } else {
      it = std::lower_bound(time.begin(), time.end(), newtime(i));
      upper = *it;
      if (it == time.begin()) {
	pos = 0;
      } else {
	pos = int(it-time.begin());
	if (type==0 && std::fabs(newtime(i)-time(pos-1)) < std::fabs(newtime(i)-time(pos))) pos -= 1;
      }
      }
      if (type==2 && newtime(i)<upper) pos--;
      idx(i) = pos;
    }
    return(idx);
  }
/*}}}*/

  arma::mat interpolate(const arma::mat &input,const double tau,const bool locf) {/*{{{*/
    vec time = input.col(0);
    unsigned n = time.n_elem;
    double t0 = time(0);
    double tn = time(n-1);
    unsigned N = std::ceil((tn-t0)/tau)+1;
    mat input2(N, input.n_cols);
    unsigned cur = 0;
    input2.row(0) = input.row(0);
    double curtime = t0;
    rowvec slope(input.n_cols);
    if (locf) {
      slope.fill(0); slope(0) = 1;
    } else {
      slope = (input.row(cur+1)-input.row(cur))/(time(cur+1)-time(cur));
    }
    for (unsigned i=0; i<N-1; i++) {
      while (time(cur+1)<curtime) {
      cur++;
      if (cur==(n-1)) break;
      if (!locf)
	slope = (input.row(cur+1)-input.row(cur))/(time(cur+1)-time(cur));
      }
      double delta = curtime-time(cur);
      input2.row(i) = input.row(cur) + slope * delta;
      curtime += tau;
    }
    // double tau2 = tn-input2(N-2,0);
    input2.row(N-1) = input.row(input.n_rows-1);
    return( input2 );
  }
/*}}}*/

arma::uvec pmini(const arma::uvec& y,const int N) {/*{{{*/
    arma::uvec res  = y ;            // residuals
    for(int i=0; i<y.n_elem; ++i) 
	    if (y(i) < N) res(i) = y(i); else res(i) = N;
    return(res); 
}
/*}}}*/

arma::colvec pminv(const arma::colvec& y,const double  N) {/*{{{*/
    arma::colvec res  = y ;            // residuals
    for(int i=0; i<y.n_elem; ++i) 
	    if (y(i) < N) res(i) = y(i); else res(i) = N;
    return(res); 
}
/*}}}*/

  arma::mat interpolate2(const arma::mat& input, const arma::colvec& tau) {/*{{{*/
    vec time = input.col(0);
    vec ftime = input.col(1);
    unsigned n = time.n_elem;
    unsigned N = tau.n_elem; 
    mat output(N, input.n_cols);

    arma::uvec pos=approx(time, tau, 2);
    pos=pmini(pos,n-1); 
    arma::uvec posp1=pmini(pos+1,n-1); 

    output.col(0) = tau; 
    for(int i=0;i<N;++i) 
    {
        if (posp1(i)==pos(i)) output(i,1) = ftime(n-1); 
        else  {
	    double slope=
	    (ftime(posp1(i))-ftime(pos(i)))/(time(posp1(i))-time(pos(i)));
             output(i,1) = ftime(pos(i)) + slope*(tau(i)-time(pos(i)));
        }
    }
    return(output);
  }
/*}}}*/

  double interpolate3(const arma::mat& input, const double dtau) {/*{{{*/
    vec time = input.col(0);
    vec ftime = input.col(1);
    unsigned n = time.n_elem;
    unsigned N = 1;
    double output; 
    vec tau(1);
    tau(0)=dtau; 

    arma::uvec pos=approx(time, tau, 2);
    pos=pmini(pos,n-1); 
    arma::uvec posp1=pmini(pos+1,n-1); 

    for(int i=0;i<N;++i) 
    {
        if (posp1(i)==pos(i)) output = ftime(n-1); 
        else  {
	    double slope=
	    (ftime(posp1(i))-ftime(pos(i)))/(time(posp1(i))-time(pos(i)));
             output = ftime(pos(i)) + slope*(tau(i)-time(pos(i)));
        }
    }
    return(output);
  }
/*}}}*/

arma::colvec  ilapC(const double theta, const arma::colvec t)
{/*{{{*/
colvec res = t;
double itheta=1/theta; 
res= (pow(t,-itheta)-1)/(itheta);
return(res);
}/*}}}*/

arma::colvec simbase(const arma::mat& cum,const double rr,const double max,const double entry,const  int maxit)
{/*{{{*/
	arma::mat cumi=cum;
	cumi.col(0)=cum.col(1); 
	cumi.col(1)=cum.col(0); 
	double lentry=entry; 

	arma::colvec sims(maxit); 
	Rcpp::NumericVector texp= Rcpp::rexp(maxit, rr); 
	int ant=0; 

	for (unsigned i=0; i<maxit; i++) {
	      double cumentry=interpolate3(cum,lentry); 
	      double t1= cumentry+texp(i);
	      double rrx=interpolate3(cumi,t1); 
	      sims(i)=rrx; 
	      lentry=rrx; 
	      ant=ant+1; 
	      if (rrx>=max) break; 
	}
	arma::colvec simss=sims.subvec(0,ant-1); 

return(simss);
}
/*}}}*/

// [[Rcpp::export(name=".rchazC")]]
arma::colvec rchazC(const arma::mat& cum,const arma::colvec rr,const arma::colvec entry)
{/*{{{*/
	arma::mat cumi=cum;
	cumi.col(0)=cum.col(1); 
	cumi.col(1)=cum.col(0); 
        unsigned N = rr.n_elem; 

//        Rcpp::NumericVector texp= Rcpp::rexp(N, 1); 
	colvec texp= Rcpp::rexp(N, 1); 
	colvec cumentry=interpolate2(cum,entry).col(1); 
	colvec t1=texp/rr+cumentry;
	colvec rrx=interpolate2(cumi,t1).col(1); 
return(rrx);
}
/*}}}*/

// [[Rcpp::export(name=".simGL")]]
arma::mat simGL(const arma::mat& dcum,const  arma::colvec& St,const  arma::colvec& rr,
		const  arma::colvec& rd, const 	arma::colvec& z,const arma::colvec& fz,
		const  arma::colvec& tc,
		const int type, const double theta,const  int maxit,const double share)
{/*{{{*/
    unsigned n = rr.n_elem;
     mat basei=dcum;
     mat sims(maxit*n,4);
     colvec dbase1=dcum.col(1); 
     colvec base1=dcum.col(1); 

     int start=0; 
	for (unsigned i=0; i<n; i++) {
               if (type==1) base1=fz[i]*cumsum(dbase1/pow(St,z[i]*rd[i]));  
               if (type==2)  {
//                  colvec Stt =exp(-z[i]*ilapC(1/theta,pow(St,rd[i])));
                  colvec Stt = pow(St,rd[i]);
                  colvec gtt =exp(theta*log(St)*rd[i]);
//                  base1=fz[i]*cumsum((dbase1/Stt)); 
                  base1=fz[i]*cumsum(dbase1/(Stt%gtt))/share; 
	       }
               if (type==3)  {
                  colvec gtt =exp(theta*log(St)*rd[i]);
                  base1=fz[i]*cumsum(dbase1/gtt)/share; 
	       }
               if (type==4)  {
                  colvec gtt =1-(1-St)*rd[i]; // model with survival on 
                  base1=fz[i]*cumsum(dbase1/gtt)/share; 
	       }
	       basei.col(1)=base1; 

	       arma::colvec simi=simbase(basei,rr(i),tc(i),0,maxit); 
	       simi=pminv(simi,tc(i)); 
//	       simi.print(); printf("==================\n"); printf(" %d \n",ni); 
	       unsigned ni=simi.n_elem; 
	       vec vi(ni,fill::value(i)); 
	       vec di(ni,fill::value(tc(i))); 
	      sims(span(start,start+ni-1),0)=vi;  
	      if (ni>=2) 
	         for (unsigned j=1; j<ni; j++) sims(start+j,1)=simi(j-1);   
              sims(span(start,start+ni-1),2)=simi;
              sims(span(start,start+ni-1),3)=di; 
              start=start+ni; 
	}
	mat outsims=sims.rows(0,start-1); 

return(outsims);
}
/*}}}*/


// [[Rcpp::export(name=".simSurvZ")]]
arma::mat simSurvZ(const arma::mat& St, const arma::colvec& rd,const arma::colvec& z,
		const double theta, const int type)
{/*{{{*/
    unsigned n = rd.n_elem;
    unsigned N = St.n_rows;
    mat basei=St;
    colvec sims(n);
    colvec St1=St.col(1); 
    colvec Stt=St1;
    double max=St(N-1,0); 

   for (unsigned i=0; i<n; i++) {
      if (type==1) Stt=pow(St1,z(i)*rd(i));  
      if (type==2) Stt=exp(-z[i]*ilapC(1/theta,pow(St1,rd(i))));
      basei.col(1)=-log(Stt); 
//    basei.print(); printf("==== %d %d %lf \n",i,n,max); 
      arma::colvec simi=simbase(basei,1,max,0,1); 
      sims(i)=simi(0); 
   }
return(sims);
}
/*}}}*/

// [[Rcpp::export(name=".tildeLambda1")]]
arma::mat tildeLambda1(const arma::colvec& dLambda1, const arma::colvec& LambdaD, 
		const arma::colvec& r1, const arma::colvec& rd, const arma::colvec& theta, const IntegerVector id)
{/*{{{*/
    unsigned n = rd.n_elem;
    unsigned N = dLambda1.n_elem;

    mat res(N,3);
    colvec resi(n); colvec Dresi(n); colvec D2resi(n);
    colvec inc(n); colvec Dinc(n); colvec D2inc(n);
    colvec rdtheta=rd%theta; 
    colvec rd2=rd%rd; 

   for (unsigned i=0; i<N; i++) {
	   if (dLambda1(i)>0.0000000000001) { 
		   resi=r1%exp((rdtheta)*LambdaD(i))*dLambda1(i); 
		   Dresi=rd%(r1*LambdaD(i))%exp(LambdaD(i)*rdtheta)*dLambda1(i); 
		   D2resi=rd2%(r1*LambdaD(i)*LambdaD(i))%exp(LambdaD(i)*rdtheta)*dLambda1(i); 
		   inc=inc+resi;
		   Dinc=Dinc+Dresi;
		   D2inc=D2inc+D2resi;
	   }
	   res(i,0)=inc(id(i)); 
	   res(i,1)=Dinc(id(i)); 
	   res(i,2)=D2inc(id(i)); 
   }
   return(res); 
}
/*}}}*/

// [[Rcpp::export(name=".tildeLambda1R")]]
arma::mat tildeLambda1R(const arma::colvec& dLambda1, const arma::colvec& LambdaD, 
		const arma::colvec& r1, const arma::colvec& rd, const arma::colvec& theta,
		const IntegerVector id, const arma::colvec& sign)
{/*{{{*/
    unsigned n = rd.n_elem;
    unsigned N = dLambda1.n_elem;

    mat res(N,3);
    colvec resi(n); colvec Dresi(n); colvec D2resi(n);
    colvec inc(n); colvec Dinc(n); colvec D2inc(n);
    colvec atrisk(n); 
    colvec rdtheta=rd%theta; colvec rd2=rd%rd; colvec rdr1=rd%r1; colvec rd2r1=rd2%r1; 

   for (unsigned i=0; i<N; i++) {
	   if (sign(i)<0) atrisk(id(i))=1; 
	   if (dLambda1(i)>0.0000000000001) { 
		   resi=r1%exp((rdtheta)*LambdaD(i))*dLambda1(i); 
		   Dresi=rdr1%exp(LambdaD(i)*rdtheta)*dLambda1(i)*LambdaD(i);
		   D2resi=rd2r1%exp(LambdaD(i)*rdtheta)*dLambda1(i)*LambdaD(i)*LambdaD(i); 
		   inc=inc+resi%atrisk;
		   Dinc=Dinc+Dresi%atrisk;
		   D2inc=D2inc+D2resi%atrisk;
	   }
	   res(i,0)=inc(id(i)); 
	   res(i,1)=Dinc(id(i)); 
	   res(i,2)=D2inc(id(i)); 
	   if (sign(i)>0) atrisk(id(i))=0; 
   }
   return(res); 
}
/*}}}*/

//
//RcppExport SEXP simSurvZ(SEXP iSt,SEXP ird, SEXP iz, SEXP itheta, SEXP itype ) {/*{{{*/
//	arma::colvec z = Rcpp::as<arma::colvec>(iz);
//	arma::colvec rd = Rcpp::as<arma::colvec>(ird);
//	arma::mat St = Rcpp::as<arma::mat>(iSt);
//	int type = Rcpp::as<int>(itype);
//	double theta = Rcpp::as<double>(itheta);
//
//	arma::mat res = simSz(St,rd,z,theta,type); 
//
//	res.print(); 
//
//	List rres;
//	rres["res"]=res;
//
////	return(Rcpp::wrap(res)); 
//}/*}}}*/
//
//RcppExport SEXP simGL(SEXP idcum,SEXP iSt,SEXP irr,
//		SEXP  ird, SEXP iz,SEXP ifz, 
//		SEXP  itc, SEXP itype, SEXP itheta,
//		SEXP imaxit)
//{/*{{{*/
//	arma::mat dcum = Rcpp::as<arma::mat>(idcum);
//	arma::mat St = Rcpp::as<arma::mat>(iSt);
//	arma::colvec z = Rcpp::as<arma::colvec>(iz);
//	arma::colvec rd = Rcpp::as<arma::colvec>(ird);
//	arma::colvec rr = Rcpp::as<arma::colvec>(irr);
//	arma::colvec fz = Rcpp::as<arma::colvec>(ifz);
//	arma::colvec tc = Rcpp::as<arma::colvec>(itc);
//	int type = Rcpp::as<int>(itype);
//	int maxit = Rcpp::as<int>(imaxit);
//	double theta = Rcpp::as<double>(itheta);
//
//	arma::mat res = simll(dcum,St,rr,rd,z,fz,tc,type,theta,maxit); 
//	res.print(); 
//
//	List rres;
//	rres["res"]=res;
//	return(rres); 
//
////	return(Rcpp::wrap(res)); 
//}/*}}}*/
//
//  arma::uvec approx(const arma::mat &time,  // Sorted time points (Ascending)
//		  const arma::mat &newtime,
//		  unsigned type=0); // (0: nearest, 1: right, 2: left)
//
//  arma::mat interpolate(const arma::mat &input, // first column is time
//		      double tau, // Time-step
//		      bool locf=false); // Last-observation-carried forward, otherwise linear interpolation

