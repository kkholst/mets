##' Simulation Helper Functions
##'
##' Internal simulation functions used for generating data from various
##' survival, competing risks, frailty, and twin/family models. These
##' functions are primarily intended for use in examples and testing.
##'
##' @section Survival and Competing Risks:
##' \describe{
##'   \item{\code{simrchaz}}{Simulate from cumulative hazard via inverse CDF.}
##'   \item{\code{simul_cifs}}{Simulate competing risks data from cumulative incidence functions.}
##'   \item{\code{simlogitSurvd}}{Simulate survival data using logistic model.}
##'   \item{\code{kumarsim}}{Simulate competing risks data (Kumaraswamy-type).}
##'   \item{\code{kumarsimRCT}}{Simulate RCT competing risks data (Kumaraswamy-type).}
##' }
##'
##' @section Clayton-Oakes and Frailty Models:
##' \describe{
##'   \item{\code{sim_ClaytonOakes_family_ace}}{Simulate family data from Clayton-Oakes ACE model.}
##'   \item{\code{sim_ClaytonOakes_twin_ace}}{Simulate twin data from Clayton-Oakes ACE model.}
##'   \item{\code{sim_Compete_twin_ace}}{Simulate twin competing risks data with ACE frailty.}
##'   \item{\code{sim_Compete_simple}}{Simulate competing risks data with shared frailty.}
##'   \item{\code{sim_Frailty_simple}}{Simulate survival data with shared frailty.}
##'   \item{\code{sim_SurvFam}}{Simulate family survival data with shared frailty.}
##' }
##'
##' @section Binomial Twin/Family Models:
##' \describe{
##'   \item{\code{sim_BinPlack}}{Simulate paired binary data (Plackett model).}
##'   \item{\code{sim_BinFam}}{Simulate binary family data.}
##'   \item{\code{sim_BinFam2}}{Simulate binary family data (alternative parameterization).}
##'   \item{\code{sim_binClaytonOakes_family_ace}}{Simulate binary family data (Clayton-Oakes ACE).}
##'   \item{\code{sim_binClaytonOakes_twin_ace}}{Simulate binary twin data (Clayton-Oakes ACE).}
##'   \item{\code{sim_binClaytonOakes_pairs}}{Simulate binary paired data (Clayton-Oakes).}
##'   \item{\code{sim_bptwin}}{Simulate from a biprobit twin model.}
##' }
##'
##' @section Nordic Twin Studies:
##' \describe{
##'   \item{\code{sim_nordictwin}}{Simulate Nordic twin registry data.}
##'   \item{\code{sim_nordic_random}}{Simulate Nordic twin data with random effects.}
##'   \item{\code{corsim_prostate_random}}{Simulate correlated prostate cancer data with random effects.}
##' }
##'
##' @name mets-simulation
##' @aliases simrchaz simul_cifs simlogitSurvd kumarsim kumarsimRCT 
##' @aliases  sim_ClaytonOakes_family_ace sim_ClaytonOakes_twin_ace
##' @aliases  sim_Compete_twin_ace sim_Compete_simple sim_Frailty_simple
##' @aliases  sim_BinPlack sim_BinFam sim_BinFam2
##' @aliases  sim_binClaytonOakes_family_ace sim_binClaytonOakes_twin_ace
##' @aliases  sim_binClaytonOakes_pairs sim_bptwin sim_SurvFam
##' @aliases  sim_nordictwin sim_nordic_random corsim_prostate_random
##' @keywords internal
NULL
