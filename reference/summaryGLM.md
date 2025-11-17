# Reporting OR (exp(coef)) from glm with binomial link and glm predictions

Reporting OR from glm with binomial link and glm predictions

## Usage

``` r
summaryGLM(object, id = NULL, fun = NULL, ...)
```

## Arguments

- object:

  glm output

- id:

  possible id for cluster corrected standard errors

- fun:

  possible function for non-standard predictions based on object

- ...:

  arguments of estimate of lava for example level=0.95

## Author

Thomas Scheike

## Examples

``` r
data(sTRACE)
sTRACE$id <- sample(1:100,nrow(sTRACE),replace=TRUE)

model <- glm(I(status==9)~sex+factor(diabetes)+age,data=sTRACE,family=binomial)
summaryGLM(model)
#> $coef
#>                   Estimate Std.Err     2.5%   97.5%   P-value
#> (Intercept)       -6.65169 0.82284 -8.26442 -5.0389 6.278e-16
#> sex                0.25832 0.22418 -0.18107  0.6977 2.492e-01
#> factor(diabetes)1  0.63305 0.30486  0.03553  1.2306 3.785e-02
#> age                0.09591 0.01099  0.07436  0.1175 2.707e-18
#> 
#> $or
#>                      Estimate         2.5%      97.5%
#> (Intercept)       0.001291843 0.0002575172 0.00648057
#> sex               1.294751030 0.8343769033 2.00914026
#> factor(diabetes)1 1.883341432 1.0361638586 3.42317957
#> age               1.100659282 1.0771942608 1.12463545
#> 
#> $fout
#> NULL
#> 
summaryGLM(model,id=sTRACE$id)
#> $coef
#>                   Estimate Std.Err     2.5%   97.5%   P-value
#> (Intercept)       -6.65169 0.77938 -8.17925 -5.1241 1.407e-17
#> sex                0.25832 0.21740 -0.16779  0.6844 2.348e-01
#> factor(diabetes)1  0.63305 0.27780  0.08856  1.1775 2.268e-02
#> age                0.09591 0.01019  0.07593  0.1159 5.032e-21
#> 
#> $or
#>                      Estimate        2.5%      97.5%
#> (Intercept)       0.001291843 0.000280412 0.00595145
#> sex               1.294751030 0.845534854 1.98262700
#> factor(diabetes)1 1.883341432 1.092604538 3.24634836
#> age               1.100659282 1.078886728 1.12287122
#> 
#> $fout
#> NULL
#> 

nd <- data.frame(sex=c(0,1),age=67,diabetes=1)
predictGLM(model,nd)
#> $coef
#>                   Estimate Std.Err     2.5%   97.5%   P-value
#> (Intercept)       -6.65169 0.82284 -8.26442 -5.0389 6.278e-16
#> sex                0.25832 0.22418 -0.18107  0.6977 2.492e-01
#> factor(diabetes)1  0.63305 0.30486  0.03553  1.2306 3.785e-02
#> age                0.09591 0.01099  0.07436  0.1175 2.707e-18
#> 
#> $pred
#>     Estimate      2.5%     97.5%
#> p1 0.6004375 0.4494731 0.7344613
#> p2 0.6605188 0.5194651 0.7778730
#> 
```
