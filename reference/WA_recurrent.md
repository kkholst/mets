# While-Alive estimands for recurrent events

Considers the ratio of means \$\$E(N(min(D,t)))/E(min(D,t))\$\$ and the
the mean of the events per time unit \$\$E(N(min(D,t))/min(D,t))\$\$
both based on IPCW etimation. RMST estimator equivalent to Kaplan-Meier
based estimator.

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

  Event formula first covariate on rhs must be a factor giving the
  treatment

- data:

  data frame

- time:

  for estimation

- cens.code:

  of censorings

- cause:

  of events

- death.code:

  of terminal events

- trans:

  possible power for mean of events per time-unit

- cens.formula:

  censoring model, default is to use strata(treatment)

- augmentR:

  covariates for model of mean ratio

- augmentC:

  covariates for censoring augmentation

- type:

  augmentation for call of binreg, when augmentC is given default is "I"
  and otherwise "II"

- marks:

  possible marks for composite outcome situation for model for counts
  with marks

- ...:

  arguments for binregATE

## References

Nonparametric estimation of the Patient Weighted While-Alive Estimand
arXiv preprint by A. Ragni, T. Martinussen, T. Scheike Mao, L. (2023).
Nonparametric inference of general while-alive estimands for recurrent
events. Biometrics, 79(3):1749–1760. Schmidli, H., Roger, J. H., and
Akacha, M. (2023). Estimands for recurrent event endpoints in the
presence of a terminal event. Statistics in Biopharmaceutical Research,
15(2):238–248.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(hfactioncpx12)

dtable(hfactioncpx12,~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 
dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,death.code=2)
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
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err  2.5% 97.5%   P-value
#> treatment0    1.572 0.09573 1.384 1.759 1.375e-60
#> treatment1    1.453 0.10315 1.251 1.656 4.376e-45
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....   0.1185  0.1407 -0.1574 0.3943     0.4
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>    Estimate Std.Err   2.5%  97.5%   P-value
#> p1   0.8457 0.05264 0.7425 0.9488 4.411e-58
#> p2   0.7555 0.05433 0.6490 0.8619 5.963e-44
#>  
#>             Estimate Std.Err     2.5%  97.5% P-value
#> [p1] - [p2]  0.09022 0.07565 -0.05805 0.2385   0.233
#> 
#>  Null Hypothesis: 
#>   [p1] - [p2] = 0 
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err   2.5%  97.5%   P-value
#> treat0   1.0725  0.1222 0.8331 1.3119 1.645e-18
#> treat1   0.7552  0.0643 0.6291 0.8812 7.508e-32
#>  
#>                     Estimate Std.Err    2.5%  97.5% P-value
#> [treat0] - [treat1]   0.3173  0.1381 0.04675 0.5879 0.02153
#> 
#>  Null Hypothesis: 
#>   [treat0] - [treat1] = 0 
```
