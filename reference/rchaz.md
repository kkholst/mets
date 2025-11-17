# Simulation of Piecewise constant hazard model (Cox).

Simulates data from piecwise constant baseline hazard that can also be
of Cox type. Censor data at highest value of the break points.

## Usage

``` r
rchaz(
  cumhazard,
  rr,
  n = NULL,
  entry = NULL,
  cum.hazard = TRUE,
  cause = 1,
  extend = FALSE
)
```

## Arguments

- cumhazard:

  cumulative hazard, or piece-constant rates for periods defined by
  first column of input.

- rr:

  relative risk for simulations, alternatively when rr=1 specify n

- n:

  number of simulation if rr not given

- entry:

  delayed entry time for simuations.

- cum.hazard:

  specifies wheter input is cumulative hazard or rates.

- cause:

  name of cause

- extend:

  to extend piecewise constant with constant rate. Default is average
  rate over time from cumulative (when TRUE), if numeric then uses given
  rate.

## Details

For a piecewise linear cumulative hazard the inverse is easy to compute
with and delayed entry x we compute \$\$\Lambda^{-1}(\Lambda(x) +
E/RR)\$\$, where RR are the relative risks and E is exponential with
mean 1. This quantity has survival function \$\$P(T \> t \| T\>x) =
exp(-RR (\Lambda(t) - \Lambda(x)))\$\$.

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
