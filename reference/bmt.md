# The Bone Marrow Transplant Data

Bone marrow transplant data with 408 rows and 5 columns.

## Format

The data has 408 rows and 5 columns.

- cause:

  a numeric vector code. Survival status. 1: dead from treatment related
  causes, 2: relapse , 0: censored.

- time:

  a numeric vector. Survival time.

- platelet:

  a numeric vector code. Plalelet 1: more than 100 x \\10^9\\ per L, 0:
  less.

- tcell:

  a numeric vector. T-cell depleted BMT 1:yes, 0:no.

- age:

  a numeric vector code. Age of patient, scaled and centered
  ((age-35)/15).

## Source

Simulated data

## References

NN

## Examples

``` r
data(bmt)
names(bmt)
#> [1] "time"     "cause"    "platelet" "age"      "tcell"   
```
