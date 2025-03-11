# mets 1.3.6
  - Development version
  - While-alive estimation : `WA_recurrent` 
  - Marks for medical cost models : `recreg` `recregIPCW`
  - New default augmentation for `binreg` `resmeanIPCW` type="II", type="I" simple outcome IPCW
  - New version of `recreg` and `cifreg` for testing and comparison with old versions `recregO` and `cifregO`
     - plot, summary, predict for these functions

# mets 1.3.5
  - `sim.phreg` and `sim.recurrent` for simulations 
  - Combining `binreg`, `binregATE` with `resmeanIPCW` and `resmeanATE`

# mets 1.3.4
  - Maintenance release
  
# mets 1.3.3
  - Lu-Tsiatis efficient logrank test and dynamic censoring augmentation: `phreg_rct`
  - Inverse Probability of treatment weighted Cox model: `phreg_IPTW`
  - Twostage randomization for survival outcome: `binregTSR` 

# mets 1.3.2
  - Extension of `recreg` (Ghosh-Lin model) to deal with composite outcomes.
  - Recurrent events regression with IPCW adjustment for fixed time point: `recregIPCW`
  - Efficient Ghosh-Lin modelling using dynamic regression augmentation.

# mets 1.3.1
  - Maintenance release

# mets 1.3.0
  - Efficient IPCW for binary data: `Effbinreg` 
  - IPCW restricted mean survival regression: `resmeanIPCW` 
  - Lin-Ying additive hazards model with fast robust standard errors: `aalenMets`  
  - mediator weighted survival mediation with robust standard errors: `mediatorSurv`
  - Examples updated
  - `dutility` function no longer casts warnings when handling formulas
  - Efficient estimation of recurrent events mean: `recurrentMarginalAIPCW` 
  - Average treatment effect for competing risks and binary data: `logitATE`, `binregATE`
  - Recurrent events regression with IPCW adjustment (Ghosh-Lin model) : `recreg`

# mets 1.2.8.1
  - Maintenance release

# mets 1.2.8
  - Augmentation of binomial regression model: `BinAugmentCifstrata`
  - Augmentation of Fine-Gray model: `FG_AugmentCifstrata`
  - Double Fine-Gray model for two causes.
  - Likelihood evaluation of `mvn` uses Moore-Penrose pseudo-inverese (threshold
    set via `lava.options(itol=...)`
  - Vignette updates

# mets 1.2.7.1
  - Maintenance release

# mets 1.2.6
  - Cumulative incidence regression `cifreg` function
    - Fine-Gray model with cloglog link of (1-F_1(t,x))
    - Logit link
  - Prototype of wildbootstrap for Cox regression with
    - confidence bands for baseline
    - with confidenence bands cumulative incidence for two cox's
  - Piecewise constant hazard: `rpch`, `ppch`
  - Test-version of multinomial regression model (via phreg): `mlogit`
  - Simulation for illness-death model: `simMultistate`
  - Haplotype modelling for discrete time-to-pregnancy models: `haplo.surv.discrete`
  - Interval censoring for discrete time logit-survival model: `interval.logitsurv.discrete`
  - Binomial Regression for competing risks data with censoring and one time point only: `binreg`

# mets 1.2.5
  - Updated predict function for `phreg`
    - with plotting functionality
    - with robust standard errors
  - New vignettes started
    - phreg robust se's for marginal Cox model
    - twostage survival  model
    - multivariate competing risks
    - recurrent events
  - `logitSurv` for fitting semiparametric proportional odds model
    - gof
    - robust standard errors for clustered case
  - twostageMLE for fast twostage fitting for clustered survival data
    with robust standard errors.
  - standard errors for twostage models now also with uncertainty from Cox baseline
  - cumulative score process test gof now also for marginal Cox models

# mets 1.2.4
  - functions `km` (Kaplan-Meier) and `cif` (cumulative incidence
    probability) with robust standard errors.
  - computation of probability of exceeding "k" events for recurrents
    processs
  - computation of probability of exceeding "k1" and "k2" events for
    bivariate recurrents processseses
  - dspline simple spline decomposition on a data frame
  - rmvn, dmvn: RNG and density for multivariate normal distribution
    with varying correlation coefficients.

# mets 1.2.3.1
  - starting values updated for twinlm method

# mets 1.2.3
  - twinlm now supports ordinal outcomes
  - optimized strata calculations in phreg
  - optimized robust standard errors in phreg
  - weights and offsets in phreg
  - weights argument added to lifetable
  - gof of phreg with fast cumulative residuals (Lin, Wei, Ying)
  - graphical gof of phreg
  - recurrent events function for marginal mean with standard errors
  - simulating recurrent events with possibly two recurrent events  and death
  - covariance calculation for recurrent events data and related bootstrap

# mets 1.2.2
  - Vignettes updated
  - Compatibility with lava version 1.5

# mets 1.2.1
  - New documentation/vignettes
  - Additional examples and unit tests
  - lifecourse plot function: lifecourse
  - block sampling function: dsample

# mets 1.2.0
  - Namespace cleaning (twostage)...
  - Dependency on R>=3.3 radix algorithm
  - Case-Control sampling for twostage model.
  - Two-stage additive gamma survival model. Additive random effects for two-stage survival
    model via pairwise composite likelihood. Simulation of family ace survival model. Function
    for computing Kendall's tau for pairs with additive gamma random effects model via simulations.
  - Two-stage additive gamma binomial model. Additive random effects for binomial model
    via pairwise composite likelihood. Simulation of family ace model.
    Function for computing pairwise concordance for for pairs with additive gamma random effects model.
  - Updated divide.conquer
  - Extra unit tests
  - force.same.cens argument with IPWC methods
  - New utility functions for data.frames
    Data processing
    - dsort
    - dreshape
    - dcut
    - drm, drename, ddrop, dkeep, dsubset
    - drelevel
    - dlag
    - dfactor, dnumeric
    Data aggregation
    - dby, dby2
    - dscalar, deval, daggregate
    - dmean, dsd, dsum, dquantile, dcor
    - dtable, dcount
    Data summaries
    - dhead, dtail,
    - dsummary,
    - dprint, dlist, dlevels, dunique

# mets 1.1.1
  - Support for left-truncation

# mets 1.1.0
  - fast.approx with 'type' argument
  - scoreMV
  - lifetable updated and new survpois function (piecewise constant
    hazard)

# mets 1.0
  - New functions biprobit.time, binomial.twostage.time.
    Automatically samples time points (approximately equidistant) up to
    last double jump time. Intial support for left truncation.
    contrast argument added to biprobit.time.
  - ipw removed (from namespace)
  - biprobit optimized for tabular data (non-continuous
    covariates). Regression design for dependence parameter (tetrachoric
    correlation) now possible.
  - predict method implemented for biprobit
  - arc-sinus transformation used for probability estimates
  - updated output of bptwin with relative recurrence risk + log-OR estimates
  - iid method for bptwin (influence function)
  - survival probabilities and start and end of intervals added to
    lifetable
  - new function 'jumptimes' for extracting jump times and possibly sample (equidistant)
  - fast.pattern updated to handle more than two categories
  - demos added to the mets package
  - divide.conquer function, folds function

# mets 0.2.8
  - Normal orthant probabilities via 'pmvn' (vectorized)
  - Parametric proportional hazards models via 'phreg.par'
  - twinlm.time function for censored twin data. Wraps the 'ipw'
    function that now also supports parametric survival models via
    phreg.par. 'grouptable' for tabulating twin-data.
  - Relative recurrence risk ratios now reported with bptwin/twinlm.
  - Grandom.cif more stable

# mets 0.2.7
  - Adapted to changes in 'timereg::comp.risk'
  - cluster.index with 'mat' argument for stacking rows of a matrix
    according to cluster-variable
  - New lava-estimator: 'normal', for ordinal data (cumulative probit)
  - fast.reshape more robust. Now also supports 'varying arguments of
    the type 'varying=-c(...)' choosing everything except '...'.

# mets 0.2.6
  - C++ source code cleanup
  - Optimization of fast.reshape

# mets 0.2.5
  - New datasets: dermalridges, dermalridgesMZ
  - Grouped analysis updated in twinlm (e.g. sex limitation model)
  - Confidence limits for genetic and environmental effects are now
    based on standard (symmetric) Wald confidence limits. (use the 'transform'
    argument of the summary method to apply logit-transform)
  - Improved output in twinlm

# mets 0.2.4
  - fast.reshape :labelnum option for both wide and long format (see
    example)
  - Compilation flags removed from Makevars files

# mets 0.2.3
  - fast.reshape bug-fix (column names)

# mets 0.2.2
  - Updated twinlm. bptwin: OS analysis
     - Better starting values for twinlm
  - Fixed claytonaokes.cpp
  - New fast cox ph regression: phreg
  - Updated two-stage estimator
  - Improved fast.reshape

# mets 0.2.1
  - fast.reshape
  - easy.binomial.twostage

# mets 0.1.4
  - Fixed cor.cpp
  - New datasets: twinstut, twinbmi, prtsim

# mets 0.1.3
  - twinlm moved to mets package, and wraps the bptwin function

# mets 0.1.2
  - code clean-up and minor bug-fixes

# mets 0.1.1
  - Random effects CIF models moved from MultiComp to mets
  - new data sets: np, multcif
  - Documentation via roxygen2
  - bug fixes

# mets 0.1.0
  - Initialization of the new package 'mets' with implementation of the
    Clayton-Oakes model with piecewise constant marginal hazards, and
    the bivariate probit random effects model (Liability model) for
    twin-data.
