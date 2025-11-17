# Simple linear spline

Constructs simple linear spline on a data frame using the formula syntax
of dutils that is adds (x-cuti)\* (x\>cuti) to the data-set for each
knot of the spline. The full spline is thus given by x and spline
variables added to the data-set.

## Usage

``` r
dspline(
  data,
  y = NULL,
  x = NULL,
  breaks = 4,
  probs = NULL,
  equi = FALSE,
  regex = mets.options()$regex,
  sep = NULL,
  na.rm = TRUE,
  labels = NULL,
  all = FALSE,
  ...
)
```

## Arguments

- data:

  if x is formula or names for data frame then data frame is needed.

- y:

  name of variable, or fomula, or names of variables on data frame.

- x:

  name of variable, or fomula, or names of variables on data frame.

- breaks:

  number of breaks, for variables or vector of break points,

- probs:

  groups defined from quantiles

- equi:

  for equi-spaced breaks

- regex:

  for regular expressions.

- sep:

  seperator for naming of cut names.

- na.rm:

  to remove NA for grouping variables.

- labels:

  to use for cut groups

- all:

  to do all variables, even when breaks are not unique

- ...:

  Optional additional arguments

## Author

Thomas Scheike

## Examples

``` r
data(TRACE)
TRACE <- dspline(TRACE,~wmi,breaks=c(1,1.3,1.7))
cca <- survival::coxph(Surv(time,status==9)~age+vf+chf+wmi,data=TRACE)
cca2 <- survival::coxph(Surv(time,status==9)~age+wmi+vf+chf+
                           wmi.spline1+wmi.spline2+wmi.spline3,data=TRACE)
anova(cca,cca2)
#> Analysis of Deviance Table
#>  Cox model: response is  Surv(time, status == 9)
#>  Model 1: ~ age + vf + chf + wmi
#>  Model 2: ~ age + wmi + vf + chf + wmi.spline1 + wmi.spline2 + wmi.spline3
#>    loglik  Chisq Df Pr(>|Chi|)    
#> 1 -6529.8                         
#> 2 -6519.6 20.323  3  0.0001455 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

nd <- data.frame(age=50,vf=0,chf=0,wmi=seq(0.4,3,by=0.01))
nd <- dspline(nd,~wmi,breaks=c(1,1.3,1.7))
pl <- predict(cca2,newdata=nd)
plot(nd$wmi,pl,type="l")

```
