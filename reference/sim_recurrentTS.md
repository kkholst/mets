# Simulate recurrent events from a two-stage model with structured gamma frailties

Simulates recurrent event data with two event types and a terminal
event, using a parametric two-stage frailty model. The construction
ensures that the marginal rates are approximately correct: conditional
on survival, \\E(dN_j \mid D \> t) \approx\\ `cumhazj`, and the hazard
of death equals `death.cumhaz`.

## Usage

``` r
sim_recurrentTS(
  n,
  cumhaz,
  cumhaz2,
  death.cumhaz = NULL,
  nu = rep(1, 3),
  share1 = 0.3,
  vargamD = 2,
  vargam12 = 0.5,
  gap.time = FALSE,
  max.recurrent = 100,
  cens = NULL,
  ...
)
```

## Arguments

- n:

  Number of subjects to simulate.

- cumhaz:

  Two-column matrix `(time, cumhaz)` giving the target marginal
  cumulative rate of the first recurrent event type.

- cumhaz2:

  Two-column matrix `(time, cumhaz)` giving the target marginal
  cumulative rate of the second recurrent event type.

- death.cumhaz:

  Two-column matrix `(time, cumhaz)` giving the cumulative hazard of the
  terminal event.

- nu:

  Numeric vector of length 3: the powers \\(\nu_1, \nu_2, \nu_3)\\
  applied to the frailty components (see Details). Must satisfy \\\nu_j
  \> -1/\text{shape}\\. Default is `rep(1, 3)`.

- share1:

  Proportion of the total death frailty variance assigned to the first
  component \\Z\_{d1}\\. The remainder goes to \\Z\_{d2}\\. Must be in
  \\(0, 1)\\. Default is `0.3`.

- vargamD:

  Total variance of the death frailty \\Z\_\text{death}\\. Default is
  `2`.

- vargam12:

  Variance of the shared recurrent-event frailty \\Z\_{12}\\. Default is
  `0.5`.

- gap.time:

  Logical. If `TRUE`, event times are drawn as gap times rather than
  calendar times. Default is `FALSE`.

- max.recurrent:

  Maximum number of recurrent events per subject. Default is `100`.

- cens:

  Rate of exponential censoring. If `NULL` (default), no administrative
  censoring is applied.

- ...:

  Further arguments passed to lower-level simulation functions.

## Value

A data frame in counting-process format (one row per event interval per
subject) with columns `id`, `start`, `stop`, `entry`, `time`, `status`,
and `death`. Attributes store the (possibly adjusted) cumulative hazards
used in simulation and the frailty parameters.

## Details

The frailty structure uses three gamma random variables \\Z\_{d1}\\,
\\Z\_{d2}\\, \\Z\_{12}\\ to induce dependence: \$\$Z\_\text{death} =
Z\_{d1} + Z\_{d2}, \quad Z_1 = Z\_{d1}^{\nu_1} Z\_{12}, \quad Z_2 =
Z\_{d2}^{\nu_2} Z\_{12}^{\nu_3}.\$\$ The parameters `share1` and
`vargamD` control how the death frailty splits between the two
components, and `vargam12` controls the shared recurrent-event frailty.
Setting \\\nu = (1,1,1)\\ with `share1 = 0.5` gives a symmetric
structure; varying \\\nu\\ allows asymmetric dependence.

## See also

[`sim_recurrentII`](http://kkholst.github.io/mets/reference/sim_recurrentII.md),
[`sim_recurrent_ts`](http://kkholst.github.io/mets/reference/sim_recurrent_ts.md)

## Author

Thomas Scheike

## Examples

``` r
data(CPH_HPN_CRBSI)
dr    <- CPH_HPN_CRBSI$terminal
base1 <- CPH_HPN_CRBSI$crbsi
base4 <- CPH_HPN_CRBSI$mechanical

rr <- sim_recurrentTS(1000, base1, base4, death.cumhaz = dr)
dtable(rr, ~death + status)
#> 
#>       status    0    1    2
#> death                      
#> 0             138 2394  384
#> 1             861    0    0
mets:::showfitsim(causes = 2, rr, dr, base1, base4)


```
