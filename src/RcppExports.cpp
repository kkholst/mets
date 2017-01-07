// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/mets.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;


// validate (ensure exported C++ functions exist before calling them)
static int mets_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP mets_RcppExport_registerCCallable() { 
    R_RegisterCCallable("mets", "mets_RcppExport_validate", (DL_FUNC)mets_RcppExport_validate);
    return R_NilValue;
}
