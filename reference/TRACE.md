# The TRACE study group of myocardial infarction

The TRACE data frame contains 1877 patients and is a subset of a data
set consisting of approximately 6000 patients. It contains data relating
survival of patients after myocardial infarction to various risk
factors.

## Format

This data frame contains the following columns:

- id:

  a numeric vector. Patient code.

- status:

  a numeric vector code. Survival status. 9: dead from myocardial
  infarction, 0: alive, 7: dead from other causes.

- time:

  a numeric vector. Survival time in years.

- chf:

  a numeric vector code. Clinical heart pump failure, 1: present, 0:
  absent.

- diabetes:

  a numeric vector code. Diabetes, 1: present, 0: absent.

- vf:

  a numeric vector code. Ventricular fibrillation, 1: present, 0:
  absent.

- wmi:

  a numeric vector. Measure of heart pumping effect based on ultrasound
  measurements where 2 is normal and 0 is worst.

- sex:

  a numeric vector code. 1: female, 0: male.

- age:

  a numeric vector code. Age of patient.

## Source

The TRACE study group.

Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L., F. E.
Nielsen, Melchior, T., Joen, T. and P. K. Andersen (1997), Does
in-hospital ventricular fibrillation affect prognosis after myocardial
infarction?, European Heart Journal 18, 919–924.

## Details

sTRACE is a subsample consisting of 300 patients.

tTRACE is a subsample consisting of 1000 patients.

## Examples

``` r
data(TRACE)
names(TRACE)
#> [1] "id"       "wmi"      "status"   "chf"      "age"      "sex"      "diabetes"
#> [8] "time"     "vf"      
```
