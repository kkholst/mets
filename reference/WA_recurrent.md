# While-Alive Estimands for Recurrent Events

Computes the "While-Alive" estimands for recurrent events in the
presence of a terminal event (death). These estimands address the
challenge of defining meaningful treatment effects when death prevents
further observation of recurrent events.

## Usage

``` r
WA_recurrent(
  formula,
  data,
  time = NULL,
  cens.code = 0,
  cause = 1,
  death.code = 2,
  trans = NULL,
  cens.formula = NULL,
  augmentR = NULL,
  augmentC = NULL,
  type = NULL,
  marks = NULL,
  ...
)
```

## Arguments

- formula:

  Formula with an `Event` object. The first covariate on the RHS must be
  a factor representing the treatment group. Can include `cluster(id)`.

- data:

  Data frame containing all variables referenced in the formula.

- time:

  Time point \\t\\ for estimation. If NULL, defaults to the maximum
  event time.

- cens.code:

  Numeric code for censoring (default 0).

- cause:

  Numeric code for the recurrent event of interest (default 1).

- death.code:

  Numeric code for the terminal event/death (default 2).

- trans:

  Power transformation for the mean of events per time-unit (default
  NULL, i.e., linear).

- cens.formula:

  Formula for the censoring model. Default is `~strata(treatment)`.

- augmentR:

  Formula for covariate augmentation in the randomization model (e.g.,
  `~age+sex`). Improves efficiency.

- augmentC:

  Formula for covariate augmentation in the censoring model. Enables
  double robustness.

- type:

  Type of augmentation for the binomial regression call. Default is "I"
  if `augmentC` is given, otherwise "II".

- marks:

  Optional marks for composite outcome situations (e.g., distinguishing
  event types in a composite endpoint).

- ...:

  Additional arguments passed to `binregATE`.

## Value

An object of class `"WA"` containing:

- RAW:

  List of raw estimates: RMST, mean number of events, ratio of means,
  and their log-transformed versions with standard errors.

- ET:

  List of estimated treatment effects: risk difference for the mean rate
  (`riskDR`) and optionally the augmented version (`riskDRC`).

- time:

  The time point used for estimation.

- cause, death.code, cens.code:

  Codes used.

- augmentR, augmentC:

  Formulas used for augmentation.

The object includes influence functions (IID) for all estimators,
allowing for further variance calculations or combination with other
estimators.

## Details

The function estimates two primary quantities:

1.  **Ratio of Means**: \$\$E(N(\min(D,t))) / E(\min(D,t))\$\$ The
    expected number of events up to time \\t\\ (censored by death \\D\\)
    divided by the expected time alive up to \\t\\.

2.  **Mean of Events per Time Unit**: \$\$E(N(\min(D,t)) /
    \min(D,t))\$\$ The expected rate of events per unit of time alive.

Estimation is based on Inverse Probability of Censoring Weighting (IPCW)
to handle administrative censoring and death. The method can be
augmented with covariates (double robust estimation) to improve
efficiency and robustness.

## References

Ragni, A., Martinussen, T., & Scheike, T. H. (2023). Nonparametric
estimation of the Patient Weighted While-Alive Estimand. arXiv preprint.

Mao, L. (2023). Nonparametric inference of general while-alive estimands
for recurrent events. Biometrics, 79(3), 1749–1760.

Schmidli, H., Roger, J. H., & Akacha, M. (2023). Estimands for recurrent
event endpoints in the presence of a terminal event. Statistics in
Biopharmaceutical Research, 15(2), 238–248.

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
dtable(hfactioncpx12,~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 
dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),data=hfactioncpx12,
                   time=2,death.code=2)
summary(dd)
#> While-Alive summaries:  
#> 
#> RMST,  E(min(D,t)) 
#>            Estimate Std.Err  2.5% 97.5% P-value
#> treatment0    1.859 0.02108 1.817 1.900       0
#> treatment1    1.924 0.01502 1.894 1.953       0
#>  
#>                           Estimate Std.Err    2.5%    97.5% P-value
#> [treatment0] - [treat.... -0.06517 0.02588 -0.1159 -0.01444  0.0118
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 6.3405, df = 1, p-value = 0.0118
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err  2.5% 97.5%   P-value
#> treatment0    1.572 0.09573 1.384 1.759 1.375e-60
#> treatment1    1.453 0.10315 1.251 1.656 4.376e-45
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....   0.1185  0.1407 -0.1574 0.3943     0.4
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 0.7085, df = 1, p-value = 0.4
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>    Estimate Std.Err   2.5%  97.5%   P-value
#> p1   0.8457 0.05264 0.7425 0.9488 4.411e-58
#> p2   0.7555 0.05433 0.6490 0.8619 5.963e-44
#>  
#>             Estimate Std.Err     2.5%  97.5% P-value
#> [p1] - [p2]  0.09022 0.07565 -0.05805 0.2385   0.233
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] - [p2] = 0 
#>  
#> chisq = 1.4222, df = 1, p-value = 0.233
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err   2.5%  97.5%   P-value
#> treat0   1.0725  0.1222 0.8331 1.3119 1.645e-18
#> treat1   0.7552  0.0643 0.6291 0.8812 7.508e-32
#>  
#>                     Estimate Std.Err    2.5%  97.5% P-value
#> [treat0] - [treat1]   0.3173  0.1381 0.04675 0.5879 0.02153
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treat0] - [treat1] = 0 
#>  
#> chisq = 5.2837, df = 1, p-value = 0.02153
```
