# CALGB 8923, twostage randomization SMART design

Data from CALGB 8923

## Format

Data from smart design id: id of subject status : 1-death, 2-response
for second randomization, 0-censoring A0 : treament at first
randomization A1 : treament at second randomization At.f : treament
given at record (A0 or A1) TR : time of response sex : 0-males,
1-females consent: 1 if agrees to 2nd randomization, censored if not R:
1 if response trt1: A0 trt2: A1

## Source

https://github.com/ycchao/code_Joint_model_SMART

## Examples

``` r
data(calgb8923)
```
