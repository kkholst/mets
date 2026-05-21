# Logrank-type test for comparing recurrent event marginal means between groups

Tests whether the marginal mean number of recurrent events differs
across groups (strata), extending the classical logrank test to the
setting of recurrent events with a competing terminal event. The test
statistic is \$\$z = \int_0^\tau w(s)\bigl\[d\hat\mu_1(s) -
d\hat\mu_2(s)\bigr\],\$\$ where \\w(s)\\ is a weight function and
\\\hat\mu_j(s)\\ is the estimated marginal mean for group \\j\\.
Variance is estimated robustly via the influence functions of Ghosh and
Lin (2000).

## Usage

``` r
test_logrankRecurrent(
  recurrent,
  death,
  weight = c("I", "II"),
  km = TRUE,
  start = 0,
  stop = NULL,
  at.risk = 5,
  cluster.id = NULL,
  ...
)
```

## Arguments

- recurrent:

  Either a `"recurrent"` object returned by
  [`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md),
  or a `"phreg"` object for the recurrent event model (in which case
  `death` must also be supplied).

- death:

  A `"phreg"` object for the terminal event model. Required when
  `recurrent` is a `"phreg"` object; ignored otherwise.

- weight:

  Character string specifying the weight scheme: `"I"`, `"II"`, or
  `"III"`. Default is `"I"`.

- km:

  Logical. If `TRUE` (default), the Kaplan-Meier estimator is used for
  the survival probability \\S(t)\\; otherwise the Nelson-Aalen
  estimator is used.

- start:

  Left truncation time for the integration. Default is `0`.

- stop:

  Right truncation time for the integration. Defaults to the last
  observed jump time.

- at.risk:

  Minimum combined risk-set size below which the weight is set to zero.
  Default is `5`.

- cluster.id:

  Optional vector of cluster identifiers for aggregating influence
  functions across clusters before forming the test statistic.

- ...:

  Currently unused.

## Value

An object of class `"estimate"` (from the lava package) with the
following components:

- coef:

  The weighted difference in marginal means between groups.

- se:

  Robust standard error of the test statistic.

- lower, upper:

  95% confidence interval bounds.

- p.value:

  Two-sided p-value for the null hypothesis of no difference.

The object also carries an `iid` attribute containing the subject-level
influence function decomposition of the test statistic, which can be
used for further inference or combination with other estimators.

## Details

Three weight schemes are available:

- `"I"`:

  (Default) \\w(t) = R_1(t) R_2(t) / (R_1(t) + R_2(t))\\, where \\R_j(t)
  = Y_j(t) / \hat S_j(t-)\\. Analogous to the standard logrank weight.

- `"II"`:

  \\w(t) = Y_j(t)\\, the raw risk-set size. Equivalent to using observed
  counts without survival adjustment.

- `"III"`:

  A modified weight incorporating the cumulative incidence, analogous to
  Gray's test for competing risks.

## References

Ghosh, D. and Lin, D. Y. (2000). Nonparametric analysis of recurrent
events and death. *Biometrics*, 56, 554–562.

## See also

[`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md),
`logrankRecurrentBase`

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf <- hfactioncpx12

## Test using two separate phreg models
xr <- phreg(Surv(entry, time, status == 1) ~ strata(treatment) + cluster(id), data = hf)
dr <- phreg(Surv(entry, time, status == 2) ~ strata(treatment) + cluster(id), data = hf)
out <- test_logrankRecurrent(xr, dr, stop = 5)
summary(out)
#> Call: estimate.default(f = FALSE, contrast = contrast, vcov = vcov(object), 
#>     coef = p)
#> ────────────────────────────────────────────────────────────
#>    Estimate Std.Err   2.5% 97.5% P-value
#> p1    37.98   27.04 -15.01 90.97  0.1601
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] = 0 
#>  
#> chisq = 1.9737, df = 1, p-value = 0.1601

## Equivalently, using a recurrent_marginal object directly
outN <- recurrent_marginal(Event(entry, time, status) ~ strata(treatment) + cluster(id),
                           data = hf, cause = 1, death.code = 2)
test_logrankRecurrent(outN)
#>    Estimate Std.Err   2.5% 97.5% P-value
#> p1    37.98   27.04 -15.01 90.97  0.1601
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] = 0 
#>  
#> chisq = 1.9737, df = 1, p-value = 0.1601
```
