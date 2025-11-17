# The Diabetic Retinopathy Data

The data was colleceted to test a laser treatment for delaying blindness
in patients with dibetic retinopathy. The subset of 197 patiens given in
Huster et al. (1989) is used.

## Format

This data frame contains the following columns:

- id:

  a numeric vector. Patient code.

- agedx:

  a numeric vector. Age of patient at diagnosis.

- time:

  a numeric vector. Survival time: time to blindness or censoring.

- status:

  a numeric vector code. Survival status. 1: blindness, 0: censored.

- trteye:

  a numeric vector code. Random eye selected for treatment. 1: left eye
  2: right eye.

- treat:

  a numeric vector. 1: treatment 0: untreated.

- adult:

  a numeric vector code. 1: younger than 20, 2: older than 20.

## Source

Huster W.J. and Brookmeyer, R. and Self. S. (1989) MOdelling paired
survival data with covariates, Biometrics 45, 145-56.

## Examples

``` r
data(diabetes)
names(diabetes)
#> [1] "id"     "time"   "status" "trteye" "treat"  "adult"  "agedx" 
```
