#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP Bhat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BhatAddGam(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BhatAddGamCC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP biprobit0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP biprobit2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bvncdf(SEXP, SEXP, SEXP);
extern SEXP claytonoakes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP claytonoakesbinRV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP claytonoakesR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP clusterindexdata(SEXP, SEXP, SEXP, SEXP);
extern SEXP clusterindexM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP familypairindex(SEXP, SEXP, SEXP);
extern SEXP FastApprox(SEXP, SEXP, SEXP, SEXP);
extern SEXP FastCluster(SEXP);
//extern SEXP FastCoxPL(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FastCoxPLstrata(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP FastCoxPLstrataPO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP FastCoxPLstrataAddGam(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,
		SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP );
//extern SEXP FastCoxPrep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FastCoxPrepStrata(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP FastLong2(SEXP, SEXP, SEXP, SEXP);
extern SEXP FastPattern(SEXP, SEXP, SEXP);
extern SEXP MatxCube(SEXP, SEXP, SEXP);
extern SEXP _mets_ApplyBy(SEXP, SEXP, SEXP);
extern SEXP _mets_ApplyBy2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mets_loglikMVN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mets_RcppExport_registerCCallable(void);
extern SEXP pBhat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pmvn0(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Dpmvn(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP survivalRV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RsurvivalRVCmarg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP survivalRV2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP survivalloglikeRVpairs(SEXP, SEXP, SEXP, SEXP, SEXP,
				   SEXP, SEXP, SEXP, SEXP, SEXP,
				   SEXP, SEXP, SEXP, SEXP, SEXP,
				   SEXP, SEXP, SEXP, SEXP);
extern SEXP twostageloglikebin(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP);
extern SEXP twostageloglikebinpairs(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP);
extern SEXP twostageloglikeRV(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP);
extern SEXP twostageloglikeRVpairs(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP);
extern SEXP Uhat(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP uniprobit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CubeVec(SEXP, SEXP,SEXP);
extern SEXP vecMatMat(SEXP, SEXP);
extern SEXP vecCPMat(SEXP);
extern SEXP OutCov(SEXP, SEXP);
extern SEXP MatxCube(SEXP, SEXP, SEXP);
extern SEXP Matdoubleindex(SEXP, SEXP, SEXP,SEXP,SEXP,SEXP);
extern SEXP CubeMat(SEXP, SEXP);
extern SEXP CubeMattime(SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP PropTestCox(SEXP, SEXP,SEXP,SEXP);
extern SEXP PropTestCoxClust(SEXP, SEXP, SEXP, SEXP,SEXP,SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP, SEXP,SEXP);
extern SEXP ModelMatrixTestCox(SEXP,SEXP,SEXP,SEXP,SEXP);
//extern SEXP simBandCumHazCox(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP revcumsumR(SEXP);
extern SEXP revcumsumstrataR(SEXP,SEXP, SEXP);
extern SEXP revcumsum2strataR(SEXP,SEXP, SEXP,SEXP,SEXP);
extern SEXP cumsum2strataR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP revcumsum2stratafdNR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP revcumsum2stratafdNR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP S0_FG_GcR(SEXP,SEXP,SEXP,SEXP, SEXP,SEXP,SEXP,SEXP, SEXP);
extern SEXP vecAllStrataR(SEXP,SEXP, SEXP);
extern SEXP vecAllStrataRevcumsumR(SEXP,SEXP, SEXP);
extern SEXP tailstrataR(SEXP,SEXP, SEXP);
extern SEXP headstrataR(SEXP,SEXP, SEXP);
extern SEXP indexstrataR(SEXP, SEXP,SEXP, SEXP,SEXP);
extern SEXP revcumsumstratasumR(SEXP,SEXP, SEXP);
extern SEXP revcumsumidstratasumR(SEXP,SEXP, SEXP,SEXP, SEXP);
extern SEXP revcumsumidstratasumCovR(SEXP,SEXP,SEXP, SEXP,SEXP, SEXP);
extern SEXP cumsumstrataR(SEXP,SEXP, SEXP);
extern SEXP diffstrataR(SEXP,SEXP, SEXP);
//extern SEXP RevCumsumstrataTwoR(SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP cumsumstrataPOR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP cumsumstrataDFGR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP cumsumstrataDFGRestrictR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP DLambetaR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP DLambetaDFGR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP riskstrataR(SEXP,SEXP,SEXP);
extern SEXP meanriskR(SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP wherestrataR(SEXP,SEXP,SEXP,SEXP);
extern SEXP maxminidR(SEXP,SEXP,SEXP);
extern SEXP cumsumstratasumR(SEXP,SEXP, SEXP);
extern SEXP cumsumidstratasumR(SEXP,SEXP, SEXP,SEXP, SEXP);
extern SEXP cumsumASR(SEXP,SEXP, SEXP);
extern SEXP cumsumidstratasumCovR(SEXP, SEXP,SEXP, SEXP,SEXP, SEXP);
extern SEXP covrfR(                SEXP,SEXP,SEXP, SEXP);
extern SEXP covrfstrataR(             SEXP,SEXP,SEXP,SEXP,SEXP, SEXP);
extern SEXP covrfstrataCovR(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, SEXP);
extern SEXP sumstrataR(SEXP,SEXP, SEXP);
extern SEXP XBmindex(SEXP,SEXP, SEXP);
//extern SEXP RVpairs2DIMRV(SEXP, SEXP);
//extern SEXP backfitEaEt(SEXP,SEXP, SEXP,SEXP,SEXP, SEXP, SEXP, SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Bhat",                              (DL_FUNC) &Bhat,                               6},
    {"BhatAddGam",                        (DL_FUNC) &BhatAddGam,                        14},
    {"BhatAddGamCC",                      (DL_FUNC) &BhatAddGamCC,                      17},
//  {"backfitEaEt",                       (DL_FUNC) &backfitEaEtt,                       9},
    {"biprobit0",                         (DL_FUNC) &biprobit0,                          8},
    {"biprobit2",                         (DL_FUNC) &biprobit2,                         10},
    {"bvncdf",                            (DL_FUNC) &bvncdf,                             3},
    {"claytonoakes",                      (DL_FUNC) &claytonoakes,                       9},
    {"claytonoakesbinRV",                 (DL_FUNC) &claytonoakesbinRV,                 10},
    {"claytonoakesR",                     (DL_FUNC) &claytonoakesR,                      6},
    {"clusterindexdata",                  (DL_FUNC) &clusterindexdata,                   4},
    {"clusterindexM",                     (DL_FUNC) &clusterindexM,                      5},
    {"cor",                               (DL_FUNC) &cor,                               40},
    {"CubeVec",                           (DL_FUNC) &CubeVec,                            3},
    {"CubeMat",                           (DL_FUNC) &CubeMat,                            2},
    {"CubeMattime",                       (DL_FUNC) &CubeMattime,                        9},
    {"familypairindex",                   (DL_FUNC) &familypairindex,                    3},
    {"FastApprox",                        (DL_FUNC) &FastApprox,                         4},
    {"FastCluster",                       (DL_FUNC) &FastCluster,                        1},
//    {"FastCoxPL",                         (DL_FUNC) &FastCoxPL,                          5},
    {"FastCoxPLstrata",                   (DL_FUNC) &FastCoxPLstrata,                   11},
    {"FastCoxPLstrataPO",                 (DL_FUNC) &FastCoxPLstrataPO,                 11},
    {"FastCoxPLstrataAddGam",             (DL_FUNC) &FastCoxPLstrataAddGam,             18},
//  {"FastCoxPrep",                       (DL_FUNC) &FastCoxPrep,                        6},
    {"FastCoxPrepStrata",                 (DL_FUNC) &FastCoxPrepStrata,                 11},
    {"FastLong2",                         (DL_FUNC) &FastLong2,                          4},
    {"FastPattern",                       (DL_FUNC) &FastPattern,                        3},
    {"indexstrataR",                      (DL_FUNC) &indexstrataR,                       5},
    {"MatxCube",                          (DL_FUNC) &MatxCube,                           3},
    {"meanriskR",                         (DL_FUNC) &meanriskR,                          5},
    {"wherestrataR",                      (DL_FUNC) &wherestrataR,                       4},
    {"maxminidR",                         (DL_FUNC) &maxminidR,                          3},
    {"Matdoubleindex",                    (DL_FUNC) &Matdoubleindex,                     6},
    {"_mets_ApplyBy",                     (DL_FUNC) &_mets_ApplyBy,                      3},
    {"_mets_ApplyBy2",                    (DL_FUNC) &_mets_ApplyBy2,                     8},
    {"_mets_loglikMVN",                   (DL_FUNC) &_mets_loglikMVN,                   14},
    {"_mets_RcppExport_registerCCallable",(DL_FUNC) &_mets_RcppExport_registerCCallable, 0},
    {"pBhat",                             (DL_FUNC) &pBhat,                              6},
    {"PropTestCox",                       (DL_FUNC) &PropTestCox,                        4},
    {"PropTestCoxClust",                  (DL_FUNC) &PropTestCoxClust,                  14},
    {"ModelMatrixTestCox",                (DL_FUNC) &ModelMatrixTestCox,                 5},
    {"pmvn0",                             (DL_FUNC) &pmvn0,                              5},
    {"revcumsumR",                        (DL_FUNC) &revcumsumR,                         1},
    {"revcumsumstrataR",                  (DL_FUNC) &revcumsumstrataR,                   3},
    {"revcumsum2strataR",                 (DL_FUNC) &revcumsum2strataR,                  5},
    {"revcumsum2stratafdNR",              (DL_FUNC) &revcumsum2stratafdNR,               7},
    {"cumsum2strataR",                    (DL_FUNC) &cumsum2strataR,                     7},
    {"S0_FG_GcR",                         (DL_FUNC) &S0_FG_GcR,                          9},
    {"vecAllStrataR",                     (DL_FUNC) &vecAllStrataR,                      3},
    {"vecAllStrataRevcumsumR",            (DL_FUNC) &vecAllStrataRevcumsumR,             3},
    {"riskstrataR",                       (DL_FUNC) &riskstrataR,                        3},
    {"revcumsumstratasumR",               (DL_FUNC) &revcumsumstratasumR,                3},
    {"cumsumstratasumR",                  (DL_FUNC) &cumsumstratasumR,                   3},
    {"cumsumstrataPOR",                   (DL_FUNC) &cumsumstrataPOR,                    6},
    {"cumsumstrataDFGR",                  (DL_FUNC) &cumsumstrataDFGR,                   6},
    {"cumsumstrataDFGRestrictR",          (DL_FUNC) &cumsumstrataDFGRestrictR,           7},
    {"DLambetaR",                         (DL_FUNC) &DLambetaR,                          8},
    {"DLambetaDFGR",                      (DL_FUNC) &DLambetaDFGR,                       8},
    {"cumsumidstratasumR",                (DL_FUNC) &cumsumidstratasumR,                 5},
    {"cumsumASR",                         (DL_FUNC) &cumsumASR,                          3},
    {"tailstrataR",                       (DL_FUNC) &tailstrataR,                        3},
    {"headstrataR",                       (DL_FUNC) &headstrataR,                        3},
    {"cumsumidstratasumCovR",             (DL_FUNC) &cumsumidstratasumCovR,              6},
    {"revcumsumidstratasumR",             (DL_FUNC) &revcumsumidstratasumR,              5},
    {"revcumsumidstratasumCovR",          (DL_FUNC) &revcumsumidstratasumCovR,           6},
    {"covrfR",                            (DL_FUNC) &covrfR,                             4},
    {"covrfstrataR",                      (DL_FUNC) &covrfstrataR,                       6},
    {"covrfstrataCovR",                   (DL_FUNC) &covrfstrataCovR,                    8},
//    {"RevCumsumstrataTwoR",             (DL_FUNC) &RevCumsumstrataTwoR,              5},
    {"cumsumstrataR",                     (DL_FUNC) &cumsumstrataR,                      3},
    {"diffstrataR",                       (DL_FUNC) &diffstrataR,                        3},
    {"sumstrataR",                        (DL_FUNC) &sumstrataR,                         3},
    {"Dpmvn",                             (DL_FUNC) &Dpmvn,                              5},
//    {"simBandCumHazCox",                  (DL_FUNC) &simBandCumHazCox,                 5},
    {"RsurvivalRVCmarg",                  (DL_FUNC) &RsurvivalRVCmarg,                   7},
    {"survivalRV",                        (DL_FUNC) &survivalRV,                        10},
    {"survivalRV2",                       (DL_FUNC) &survivalRV2,                       10},
    {"survivalloglikeRVpairs",            (DL_FUNC) &survivalloglikeRVpairs,            19},
    {"twostageloglikebin",                (DL_FUNC) &twostageloglikebin,                21},
    {"twostageloglikebinpairs",           (DL_FUNC) &twostageloglikebinpairs,           22},
    {"twostageloglikeRV",                 (DL_FUNC) &twostageloglikeRV,                 19},
    {"twostageloglikeRVpairs",            (DL_FUNC) &twostageloglikeRVpairs,            19},
    {"Uhat",                              (DL_FUNC) &Uhat,                               5},
    {"uniprobit",                         (DL_FUNC) &uniprobit,                          8},
    {"vecMatMat",                         (DL_FUNC) &vecMatMat,                          2},
    {"vecCPMat",                          (DL_FUNC) &vecCPMat,                           1},
    {"OutCov",                            (DL_FUNC) &OutCov,                             2},
    {"XBmindex",                          (DL_FUNC) &XBmindex,                           3},
//    {"RVpairs2DIMRV",                     (DL_FUNC) &RVpairs2DIMRV,                      2},
    {NULL, NULL, 0}
};

void R_init_mets(DllInfo *dll)
{
    R_registerRoutines(dll,
                       NULL,	       /* slot for .C */
                       CallEntries,    /* slot for .Call */
                       NULL,              /* slot for .Fortran */
                       NULL);   	       /* slot for .External */
    R_useDynamicSymbols(dll, TRUE);    /* visibility */
}
