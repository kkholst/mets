# Influence curve components for binomial regression ATE

Extracts and scales influence function (IC) components from a fitted
binomial regression ATE object or raw IC matrix. Depending on the object
class and requested \`type\`, returns one or more components combined by
column-binding.

## Usage

``` r
# S3 method for class 'binreg'
IC(x, type = c("coef", "DR", "G"), ...)
```

## Arguments

- x:

  An object of class `ATE`, or an object containing an `iid` component.

- type:

  Character vector specifying which influence curve components to
  return. Possible values are:

  coef

  :   Influence curve for regression coefficients

  DR

  :   Doubly robust ATE influence curve

  G

  :   G-computation ATE influence curve

  If `x` inherits from `ATE` and `type` is missing, the default is
  `c("DR", "G")`.

- ...:

  Currently unused.

## Value

A numeric matrix obtained by column-binding the selected influence curve
components. Each column corresponds to one element of `type`, with
column names preserved.

## Details

For objects of class `ATE`, the function extracts and rescales the
stored influence curve components:

- `x$iid` for coefficient ICs

- `x$riskDR.iid` for doubly robust ATE ICs

- `x$riskG.iid` for G-computation ICs

Each component is multiplied by its sample size (`NROW(.)`).

For non-`ATE` objects, the function returns `NROW(x$iid) * x$iid`.

## Examples

``` r
if (FALSE) { # \dontrun{
# default ATE behavior (DR + G)
IC.binreg(fit)

# only DR component
IC.binreg(fit, "DR")

# multiple components
IC.binreg(fit, c("coef", "DR"))
} # }
```
