# Estimation of probability of more that k events for recurrent events process

Estimation of probability of more that k events for recurrent events
process where there is terminal event, based on this also estimate of
variance of recurrent events. The estimator is based on cumulative
incidence of exceeding "k" events. In contrast the probability of
exceeding k events can also be computed as a counting process integral,
and this is implemented in prob.exceedRecurrent

## Usage

``` r
prob.exceed.recurrent(
  formula,
  data,
  cause = 1,
  death.code = 2,
  cens.code = 0,
  exceed = NULL,
  marks = NULL,
  cifmets = TRUE,
  all.cifs = FALSE,
  ...
)
```

## Arguments

- formula:

  formula

- data:

  data-frame

- cause:

  of interest

- death.code:

  for status

- cens.code:

  censoring codes

- exceed:

  values (if not given then all observed values)

- marks:

  may be give for jump-times and then exceed values needs to be
  specified

- cifmets:

  if true uses cif of mets package rather than prodlim

- all.cifs:

  if true then returns list of all fitted objects in cif.exceed

- ...:

  Additional arguments to lower level funtions

## References

Scheike, Eriksson, Tribler (2019) The mean, variance and correlation for
bivariate recurrent events with a terminal event, JRSS-C

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

oo <- prob.exceed.recurrent(Event(entry,time,status)~cluster(id),
        hfactioncpx12,cause=1,death.code=2)
plot(oo)

```
