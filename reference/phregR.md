# Fast Cox PH regression and calculations done in R to make play and adjustments easy

Fast Cox PH regression with R implementation to play and adjust in R
function: FastCoxPLstrataR

## Usage

``` r
phregR(formula, data, offset = NULL, weights = NULL, ...)
```

## Arguments

- formula:

  formula with 'Surv' outcome (see `coxph`)

- data:

  data frame

- offset:

  offsets for cox model

- weights:

  weights for Cox score equations

- ...:

  Additional arguments to lower level funtions

## Details

Robust variance is default variance with the summary.

influence functions (iid) will follow numerical order of given cluster
variable so ordering after \$id will give iid in order of data-set.

## Author

Klaus K. Holst, Thomas Scheike
