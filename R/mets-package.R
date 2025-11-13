#' @keywords internal
"_PACKAGE"

##' Analysis of Multivariate Events
##'
##' Implementation of various statistical models for multivariate
##' event history data. Including multivariate cumulative incidence models,
##' and bivariate random effects probit models (Liability models)
##'
##' @name mets-package
##' @author Klaus K. Holst and Thomas Scheike
##' @useDynLib mets, .registration=TRUE
##' @import stats splines Rcpp mvtnorm
##' @importFrom lava iid estimate bootstrap compare score information twostage
##'   %++% %ni% addvar<- blockdiag cancel Col confband constrain<- constraints
##'   covariance covariance<- coxWeibull.lvm devcoords distribution<- endogenous
##'   eventTime Expand getoutcome gof intercept<- Inverse kill<- latent latent<-
##'   lava.options lvm Model multigroup parameter<- pars regression regression<-
##'   revdiag trim IC expit logit
##' @importFrom survival Surv is.Surv concordance strata cluster finegray
##' @importFrom timereg two.stage predict.two.stage
##' @importFrom utils head tail getS3method glob2rx capture.output
##' @importFrom graphics matplot lines plot polygon par points abline title
##'   matlines legend mtext layout axis barplot mosaicplot
##' @importFrom methods hasArg
##' @importFrom grDevices dev.list devAskNewPage dev.interactive
##' @aliases mets mets-package
##' @examples
##'
##' ## To appear
##'
NULL

##' @export
lava::IC

##' @export
lava::iid

##' @export
lava::twostage

##' @export
lava::estimate

##' @export
lava::gof

##' @export
survival::strata

##' @export
survival::Surv

##' @export
survival::cluster

##' np data set
##'
##' @name np
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Migraine data
##'
##' @name migr
##' @docType data
##' @keywords data
NULL

##' Dermal ridges data (families)
##'
##' Data on dermal ridge counts in left and right hand in (nuclear) families
##' @name dermalridges
##' @docType data
##' @keywords data
##' @format Data on 50 families with ridge counts in left and right
##' hand for moter, father and each child. Family id in 'family' and
##' gender and child number in 'sex' and 'child'.
##' @source Sarah B. Holt (1952). Genetics of dermal ridges: bilateral
##' asymmetry in finger ridge-counts.  Annals of Eugenics 17 (1),
##' pp.211--231. DOI: 10.1111/j.1469-1809.1952.tb02513.x
##' @examples
##' data(dermalridges)
##' fast.reshape(dermalridges,id="family",varying=c("child.left","child.right","sex"))
NULL

##' Dermal ridges data (monozygotic twins)
##'
##' Data on dermal ridge counts in left and right hand in (nuclear) families
##' @name dermalridgesMZ
##' @docType data
##' @keywords data
##' @format Data on dermal ridge counts (left and right hand) in 18
##' monozygotic twin pairs.
##' @source Sarah B. Holt (1952). Genetics of dermal ridges: bilateral
##' asymmetry in finger ridge-counts.  Annals of Eugenics 17 (1),
##' pp.211--231. DOI: 10.1111/j.1469-1809.1952.tb02513.x
##' @examples
##' data(dermalridgesMZ)
##' fast.reshape(dermalridgesMZ,id="id",varying=c("left","right"))
NULL

##' CALGB 8923, twostage randomization SMART design
##'
##' Data from CALGB 8923 
##' @name calgb8923 
##' @docType data
##' @keywords data
##' @format Data from smart design
##' id: id of subject
##' status : 1-death, 2-response for second randomization, 0-censoring
##' A0 : treament at first randomization
##' A1 : treament at second randomization
##' At.f : treament given at record (A0 or A1)
##' TR : time of response 
##' sex : 0-males, 1-females 
##' consent: 1 if agrees to 2nd randomization, censored if not
##' R: 1 if response 
##' trt1: A0
##' trt2: A1
##' @source  https://github.com/ycchao/code_Joint_model_SMART
##' @examples
##' data(calgb8923)
NULL

##' ACTG175, block randmized study from speff2trial package
##'
##' Data from speff2trial
##' @name ACTG175
##' @docType data
##' @keywords data
##' @format Randomized study 
##' @source  Hammer et al. 1996, speff2trial package.
##' @examples
##' data(ACTG175)
NULL

##' hfaction, subset of block randmized study HF-ACtion from WA package  
##'
##' Data from HF-action trial slightly modified from WA package,
##' consisting of 741 nonischemic patients with baseline
##' cardiopulmonary test duration less than or equal to 12 minutes.
##' @name hfactioncpx12
##' @docType data
##' @keywords data
##' @format Randomized study 
##' status : 1-event, 2-death, 0-censoring
##' treatment : 1/0
##' @source  WA package, Connor et al. 2009
##' @examples
##' data(hfactioncpx12)
NULL


##' Menarche data set
##'
##' @name mena
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Multivariate Cumulative Incidence Function example data set
##'
##' @name multcif
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Stutter data set
##'
##' Based on nation-wide questionnaire answers from 33,317 Danish twins
##' @format
##' tvparnr: twin-pair id
##' zyg: zygosity, MZ:=mz, DZ(same sex):=dz, DZ(opposite sex):=os
##' stutter: stutter status (yes/no)
##' age: age
##' nr: number within twin-pair
##' @name twinstut
##' @docType data
##' @keywords data
NULL

##' BMI data set
##'
##' @format
##' Self-reported BMI-values on 11,411 subjects
##'
##' tvparnr: twin id
##' bmi: BMI (m/kg^2)
##' age: Age
##' gender: (male/female)
##' zyg: zygosity, MZ:=mz, DZ(same sex):=dz, DZ(opposite sex):=os
##' @name twinbmi
##' @docType data
##' @keywords data
NULL

##' Prostate data set
##'
##' @name prt
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' For internal use
##'
##' @title For internal use
##' @name npc
##' @rdname internal
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
##' @aliases plotcr npc nonparcuminc simnordic corsim.prostate
##'  alpha2kendall alpha2spear coefmat piecewise.twostage surv.boxarea
##'  faster.reshape piecewise.data
##'  simBinPlack simBinFam simBinFam2 simSurvFam corsim.prostate.random
##'  simnordic.random simCox sim
##'  grouptable jumptimes folds
##'  ace.family.design ascertained.pairs CCbinomial.twostage
##'  coarse.clust concordanceTwinACE concordanceTwostage
##'  fast.cluster force.same.cens ilap
##'  kendall.ClaytonOakes.twin.ace kendall.normal.twin.ace
##'  make.pairwise.design make.pairwise.design.competing
##'  matplot.mets.twostage object.defined p11.binomial.twostage.RV
##'  predictPairPlack simbinClaytonOakes.family.ace
##'  simbinClaytonOakes.pairs simbinClaytonOakes.twin.ace
##'  simClaytonOakes.family.ace simClaytonOakes.twin.ace simFrailty.simple
##'  simCompete.simple simCompete.twin.ace twin.polygen.design
##'  procform procform3 procformdata drop.specials
NULL


##' Rates for HPN program for patients of Copenhagen Cohort
##'
##' @name CPH_HPN_CRBSI 
##' @format
##'  crbsi: cumulative rate of catheter related bloodstream infection in HPN patients of Copenhagen
##'  mechanical: cumulative rate of Mechanical (hole/defect) complication for catheter of HPN patients of Copenhagen
##' trombo: cumulative rate of Occlusion/Thrombosis complication for catheter of HPN patients of Copenhagen
##' terminal: rate of terminal event, patients leaving the HPN program
##' @docType data
##' @keywords data
##' @source Estimated data
NULL

##' haplo fun data 
##'
##' @name haplo 
##' @format
##' hapfreqs : haplo frequencies
##' haploX:  covariates and response for haplo survival discrete survival
##' ghaplos:  haplo-types for subjects of haploX data
##' @docType data
##' @keywords data
##' @source Estimated data
NULL

##' ttpd discrete survival data on interval form
##'
##' @name ttpd
##' @docType data
##' @keywords data
##' @source Simulated data
NULL
