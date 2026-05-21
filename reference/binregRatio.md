# Percentage of Years Lost Due to a Cause Regression

Estimates the percentage of the restricted mean time lost (RMTL) that is
attributable to a specific cause and models how covariates affect this
percentage using IPCW regression.

## Usage

``` r
binregRatio(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
  type = c("III", "II", "I"),
  offset = NULL,
  weights = NULL,
  cens.weights = NULL,
  cens.model = ~+1,
  se = TRUE,
  relative.to.causes = NULL,
  kaplan.meier = TRUE,
  cens.code = 0,
  no.opt = FALSE,
  method = "nr",
  augmentation = NULL,
  outcome = c("rmtl", "cif"),
  model = c("logit", "exp", "lin"),
  Ydirect = NULL,
  ...
)
```

## Arguments

- formula:

  Formula with an outcome (see `coxph`). The first covariate on the RHS
  is typically the treatment or group indicator. Can include
  `cluster(id)`.

- data:

  Data frame containing the variables.

- cause:

  Numeric code of the cause of interest.

- time:

  Time point \\t\\ for the analysis. Required.

- beta:

  Starting values for optimization (default NULL, uses zeros).

- type:

  Type of estimator: `"I"` (IPCW only), `"II"` (IPCW + augmentation), or
  `"III"` (IPCW + complex augmentation). Default is `"III"`.

- offset:

  Offsets for the partial likelihood.

- weights:

  Weights for the score equations.

- cens.weights:

  External censoring weights (if provided, `cens.model` is ignored).

- cens.model:

  Formula for the censoring model (default `~+1`, stratified KM). Can
  include [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) for
  stratified censoring.

- se:

  Logical; if TRUE, computes standard errors based on IPCW (default
  TRUE).

- relative.to.causes:

  If not NULL, compares the RMTL of the specified `cause` to the RMTL of
  the causes in this vector (the denominator becomes the sum of these
  causes).

- kaplan.meier:

  Logical; if TRUE, uses Kaplan-Meier for IPCW weights; if FALSE, uses
  \\\exp(-\text{cumulative hazard})\\.

- cens.code:

  Censoring code (default 0).

- no.opt:

  Logical; if TRUE, skips optimization and uses `beta` directly.

- method:

  Optimization method: `"nr"` (Newton-Raphson) or `"nlm"`.

- augmentation:

  Initial augmentation term (used for type "II" and "III").

- outcome:

  Outcome type: `"rmtl"` (years lost) or `"cif"` (cumulative incidence).

- model:

  Link function: `"logit"` (default), `"exp"`, or `"lin"`.

- Ydirect:

  Matrix with two columns (numerator, denominator) to use directly as
  the response.

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"binreg"` and `"ratio"` containing:

- coef:

  Coefficient estimates.

- se.coef:

  Standard errors.

- var:

  Variance-covariance matrix.

- iid:

  Influence function decomposition (with censoring adjustment).

- iidI:

  Influence function without censoring adjustment.

- naive.var:

  Variance assuming known censoring.

- time:

  Time point used.

- cause:

  Cause of interest.

- Causes:

  Set of causes considered in the denominator.

- Yipcw:

  IPCW-adjusted response matrix.

- coefI, varI:

  Results from the initial (type "I") fit.

- augmentation:

  Final augmentation term used.

## Details

Let the total years lost be \\Y = t - \min(T, t)\\ and the years lost
due to cause 1 be \\Y_1 = I(\epsilon=1) (t - \min(T, t))\\. The function
models the ratio: \$\$ \text{logit}\left( \frac{E(Y_1 \| X)}{E(Y \| X)}
\right) = X^T \beta \$\$

Estimation is based on a binomial regression IPCW response estimating
equation: \$\$ X \left( \Delta^{\text{ipcw}}(t) \left( Y \cdot
\text{expit}(X^T \beta) - Y_1 \right) \right) = 0 \$\$ where
\\\Delta^{\text{ipcw}}(t) = I(\min(t,T) \< C) / G_c(\min(t,T))\\ is the
IPCW adjustment.

The function supports three types of estimators:

- `"I"`: Classical outcome IPCW regression (no augmentation).

- `"II"`: Adds a censoring augmentation term \\X \int E(Z(t)\|
  T\>s)/G_c(s) d \hat M_c\\ to improve efficiency (requires an initial
  estimate of \\\beta\\).

- `"III"`: Adds a more complex augmentation term separating the
  expectations of \\Y\\ and \\Y_1\\ for further efficiency gains.

The variance is based on the squared influence functions (IID). A
"naive" variance (assuming known censoring) is also provided for
comparison.

## References

Scheike, T. & Tanaka, S. (2025). Restricted mean time lost ratio
regression: Percentage of restricted mean time lost due to specific
cause. WIP.

## See also

[`resmeanIPCW`](http://kkholst.github.io/mets/reference/resmeanIPCW.md),
[`binreg`](http://kkholst.github.io/mets/reference/binreg.md)

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$time <- bmt$time+runif(408)*0.001

rmtl30 <- rmstIPCW(Event(time,cause!=0)~platelet+tcell+age, bmt, time=30, cause=1, outcome="rmtl")
rmtl301 <- rmstIPCW(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1)
rmtl302 <- rmstIPCW(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=2)

estimate(rmtl30)
#>             Estimate Std.Err    2.5%   97.5%   P-value
#> (Intercept)   2.7657 0.05104  2.6657  2.8658 0.000e+00
#> platelet     -0.3004 0.11176 -0.5194 -0.0813 7.201e-03
#> tcell        -0.1455 0.14445 -0.4286  0.1377 3.139e-01
#> age           0.2015 0.04824  0.1069  0.2960 2.961e-05
estimate(rmtl301)
#>             Estimate Std.Err    2.5%      97.5%    P-value
#> (Intercept)   2.4449 0.07291  2.3020  2.5878254 1.618e-246
#> platelet     -0.4519 0.16739 -0.7800 -0.1238266  6.940e-03
#> tcell        -0.4937 0.25146 -0.9866 -0.0008656  4.960e-02
#> age           0.2747 0.06538  0.1466  0.4028542  2.648e-05
estimate(rmtl302)
#>             Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)  1.45849  0.1362  1.1915 1.7255 9.621e-27
#> platelet    -0.01890  0.2224 -0.4548 0.4170 9.323e-01
#> tcell        0.40759  0.2648 -0.1113 0.9265 1.237e-01
#> age          0.04576  0.1191 -0.1877 0.2792 7.008e-01

## Percentage of total RMTL due to cause 1
rmtlratioI <- rmtlRatio(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1)
summary(rmtlratioI)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.973135  0.177855  0.624545  1.321725  0.0000
#> platelet    -0.402875  0.323413 -1.036753  0.231002  0.2129
#> tcell       -0.845357  0.419214 -1.667002 -0.023712  0.0437
#> age          0.214462  0.168387 -0.115571  0.544495  0.2028
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.64623 1.86740 3.7499
#> platelet     0.66840 0.35460 1.2599
#> tcell        0.42940 0.18881 0.9766
#> age          1.23920 0.89086 1.7237
#> 
#> 

newdata <- data.frame(platelet=1, tcell=1, age=1)
pp <- predict(rmtlratioI, newdata)
pp
#>        pred        se    lower     upper
#> 1 0.4848458 0.1092229 0.270769 0.6989226

## Percentage of total cumulative incidence due to cause 1
cifratio <- binregRatio(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1, model="cif")
summary(cifratio)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>               Estimate    Std.Err       2.5%      97.5% P-value
#> (Intercept)  0.7239521  0.0355963  0.6541846  0.7937196  0.0000
#> platelet    -0.0889619  0.0730427 -0.2321228  0.0541991  0.2232
#> tcell       -0.1951701  0.1005726 -0.3922888  0.0019486  0.0523
#> age          0.0452930  0.0358278 -0.0249281  0.1155142  0.2062
#> 
#> 
#> 
pp <- predict(cifratio, newdata)
pp
#>        pred        se     lower     upper
#> 1 0.4851132 0.1050915 0.2791377 0.6910887
```
