# Weibull-Cox regression

Fits a Cox-Weibull with cumulative hazard given by \$\$ \Lambda(t) =
\lambda \cdot t^s \$\$ where \\s\\ is the shape parameter, and
\\\lambda\\ the rate parameter. We here allow a regression model for
both parameters \$\$\lambda := \exp(\beta^\top X)\$\$ \$\$s :=
\exp(\gamma^\top Z)\$\$ as defined by \`formula\` and \`shape.formula\`
respectively.

## Usage

``` r
phreg_weibull(formula, shape.formula = ~1, data, control = list())
```

## Arguments

- formula:

  Formula for proportional hazards. The right-handside must be an
  \[Event\] or \[Surv\] object (with right-censoring and possibly
  delayed entry).

- shape.formula:

  Formula for shape parameter

- data:

  data.frame

- control:

  control arguments to optimization routine \[stats::nlmbin\]

## Value

\`phreg.par\` object

## Details

The parametrization

## See also

\[mets::phreg()\]

## Author

Klaus Kähler Holst, Thomas Scheike

## Examples

``` r
data(sTRACE, package="mets")
sTRACE$entry <- 0
fit1 <- phreg_weibull(Event(entry, time, status == 9) ~ age,
             shape.formula = ~age, data = sTRACE)
tt <- seq(0,10, length.out=100)
pr1 <- predict(fit1, newdata = sTRACE[1, ], times = tt)
fit2 <- phreg(Event(time, status == 9) ~ age, data = sTRACE)
pr2 <- predict(fit2, newdata = sTRACE[1, ], se = FALSE)
if (interactive()) {
   plot(pr2$times, pr2$surv, type="s")
   lines(tt, pr1[,1,1], col="red", lwd=2)
}
```
