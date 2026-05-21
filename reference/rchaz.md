# Simulation of Piecewise Constant Hazard Model (Cox)

Simulates data from a piecewise constant baseline hazard that can also
be of Cox type. Censors data at the highest value of the break points.

## Usage

``` r
rchaz(cumhazard, rr, n = NULL, entry = NULL, cause = 1, extend = FALSE)
```

## Arguments

- cumhazard:

  Cumulative hazard matrix (columns: time, cumulative hazard), or
  piece-constant rates for periods defined by the first column of input.

- rr:

  Relative risk vector for simulations, or alternatively when `rr=1`
  specify `n`.

- n:

  Number of simulations if `rr` not given.

- entry:

  Delayed entry time for simulations (optional).

- cause:

  Name/code of the event cause (default 1).

- extend:

  To extend piecewise constant with constant rate beyond the last break
  point. Default is `FALSE`. If `TRUE`, extends with average rate over
  time from cumulative. If numeric, uses the given rate.

## Value

A data frame containing:

- entry:

  Entry times.

- time:

  Event/censoring times.

- status:

  Event status (1=event, 0=censored).

- rr:

  Relative risks used.

Attributes include:

- cumhaz:

  The cumulative hazard used.

- extend.rate:

  The extension rate if used.

## Details

For a piecewise linear cumulative hazard, the inverse is easy to
compute. With delayed entry \\x\\, we compute:
\$\$\Lambda^{-1}(\Lambda(x) + E/RR)\$\$ where \\RR\\ are the relative
risks and \\E\\ is exponential with mean 1. This quantity has survival
function: \$\$P(T \> t \| T\>x) = \exp(-RR (\Lambda(t) -
\Lambda(x)))\$\$

## Author

Thomas Scheike

## Examples

``` r
chaz <-  c(0,1,1.5,2,2.1)
breaks <- c(0,10,   20,  30,   40)
cumhaz <- cbind(breaks,chaz)
n <- 10
X <- rbinom(n,1,0.5)
beta <- 0.2
rrcox <- exp(X * beta)

pctime <- rchaz(cumhaz,n=10)
pctimecox <- rchaz(cumhaz,rrcox,entry=runif(n))
```
