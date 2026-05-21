# Estimate the probability of exceeding k recurrent events by time t

Estimates \\P(N(t) \geq k)\\ as a function of time \\t\\, for a range of
thresholds \\k\\, in the presence of a terminal event (death). The
estimator is based on the cumulative incidence of "reaching \\k\\
events", treating death as a competing risk. Confidence intervals are
computed on the log or plain scale.

## Usage

``` r
prob_exceed_recurrent(
  formula,
  data,
  cause = 1,
  death.code = 2,
  cens.code = 0,
  exceed = NULL,
  marks = NULL,
  all.cifs = FALSE,
  return.data = FALSE,
  conf.type = c("log", "plain"),
  level = 0.95,
  ...
)
```

## Arguments

- formula:

  A formula with an
  [`Event`](http://kkholst.github.io/mets/reference/Event.md) response
  giving the exit time and status (and optionally entry time). The
  right-hand side may include
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html) and
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html).

- data:

  A data frame containing all variables in `formula`.

- cause:

  Integer code for the recurrent event of interest. Default is `1`.

- death.code:

  Integer code for the terminal event. Default is `2`.

- cens.code:

  Integer code for censoring. Default is `0`.

- exceed:

  Integer vector of thresholds \\k\\ to evaluate. If `NULL` (default),
  all observed cumulative counts are used.

- marks:

  Optional numeric vector of event weights. If non-`NULL`, cumulative
  counts are formed as weighted sums of events rather than simple
  counts. Must have the same length as `nrow(data)`.

- all.cifs:

  Logical. If `TRUE`, the fitted `cif` object for each threshold is
  returned in `cif.exceed`. Default is `FALSE`.

- return.data:

  Logical. If `TRUE`, the constructed dataset for each threshold is
  returned in `dataList`. Default is `FALSE`.

- conf.type:

  Type of confidence interval transformation: `"log"` (default) or
  `"plain"`.

- level:

  Confidence level. Default is `0.95`.

- ...:

  Further arguments passed to
  [`cif`](http://kkholst.github.io/mets/reference/cif.md).

## Value

An object of class `"exceed"` with the following components:

- time:

  Vector of evaluation time points.

- prob:

  Array of dimension `(length(time), length(exceed) + 1, nstrata)`
  containing \\P(N(t) \geq k)\\ for each threshold and stratum. The
  first column gives \\P(N(t) \< \text{exceed}\[1\])\\.

- se.prob:

  Standard errors of `prob`.

- lower, upper:

  Pointwise confidence interval bounds.

- meanN:

  Estimated mean number of events \\E(N(t))\\ (single stratum only;
  `NULL` for stratified analyses).

- meanN2, varN:

  Second moment and variance of \\N(t)\\ (single stratum only).

- exceed:

  Thresholds evaluated (excluding zero).

- cif.exceed:

  List of fitted `cif` objects (if `all.cifs = TRUE`).

- dataList:

  List of datasets for each threshold (if `return.data = TRUE`).

- nstrata, strata.levels, strata.name:

  Stratification information.

Use [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods for
visualisation and tabulation.

## Details

For each threshold \\k\\ in `exceed`, the function identifies the first
time each subject reaches \\k\\ events, then fits a competing risks
model ([`cif`](http://kkholst.github.io/mets/reference/cif.md)) with
"reaching \\k\\ events" as the event of interest and death as the
competing event. Strata are supported. When `marks` is `NULL`, each
event contributes equally; otherwise events are weighted by their mark
values before cumulative counts are formed.

## References

Scheike, T. H., Eriksson, L., and Tribler, P. (2019). The mean, variance
and correlation for bivariate recurrent events with a terminal event.
*Journal of the Royal Statistical Society, Series C*, 68(5).

## See also

[`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md),
[`cif`](http://kkholst.github.io/mets/reference/cif.md)

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
dtable(hfactioncpx12, ~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 

oo <- prob_exceed_recurrent(Event(entry, time, status) ~ cluster(id),
                             hfactioncpx12, cause = 1, death.code = 2)
plot(oo)

summary(oo, times = c(1, 2, 5))
#> $prob
#>      times                 N<1 exceed>=1 exceed>=2 exceed>=3  exceed>=4
#> [1,]     1 0.9978807 0.5747460 0.4252540 0.2008652 0.1012955 0.04794006
#> [2,]     2 1.9967128 0.3925156 0.6074844 0.3509483 0.2205076 0.13989818
#> [3,]     5 3.9793816 0.1925999 0.8074001 0.5477499 0.3899373 0.29900312
#>       exceed>=5  exceed>=6   exceed>=7
#> [1,] 0.03153223 0.01371012 0.008229099
#> [2,] 0.10092792 0.05533511 0.035595440
#> [3,] 0.19615192 0.14357991 0.103037717
#> 
#> $se
#>      times                  N<1  exceed>=1  exceed>=2  exceed>=3   exceed>=4
#> [1,]     1 0.9978807 0.01827977 0.01827977 0.01481729 0.01116527 0.007907269
#> [2,]     2 1.9967128 0.01862412 0.01862412 0.01832598 0.01592337 0.013395849
#> [3,]     5 3.9793816 0.02129779 0.02129779 0.02413438 0.02515176 0.024909447
#>        exceed>=5   exceed>=6   exceed>=7
#> [1,] 0.006470709 0.004305759 0.003345689
#> [2,] 0.011774638 0.008883381 0.007151615
#> [3,] 0.020706361 0.019083117 0.016923719
#> 
#> $lower
#>      times                                                             
#> [1,]     1 0.9978807 0.6091060 0.3908940 0.1738256 0.08161439 0.0346977
#> [2,]     2 1.9967128 0.4279431 0.5720569 0.3168070 0.19140637 0.1159594
#> [3,]     5 3.9793816 0.2332821 0.7667179 0.5024323 0.34362952 0.2539590
#>                                        
#> [1,] 0.02109017 0.007408243 0.003709206
#> [2,] 0.08029840 0.040397144 0.024009120
#> [3,] 0.15949140 0.110652445 0.074677237
#> 
#> $upper
#>      times                                                             
#> [1,]     1 0.9978807 0.5373658 0.4626342 0.2321109 0.1257227 0.06623635
#> [2,]     2 1.9967128 0.3548940 0.6451060 0.3887690 0.2540334 0.16877894
#> [3,]     5 3.9793816 0.1497591 0.8502409 0.5971550 0.4424855 0.35203666
#>                                     
#> [1,] 0.0471443 0.02537272 0.01825676
#> [2,] 0.1268574 0.07579680 0.05277309
#> [3,] 0.2412392 0.18630579 0.14216877
#> 
```
