# Average Treatment effect for Restricted Mean for censored competing risks data using IPCW

Under the standard causal assumptions we can estimate the average
treatment effect E(Y(1) - Y(0)). We need Consistency, ignorability (
Y(1), Y(0) indep A given X), and positivity.

## Usage

``` r
resmeanATE(formula, data, model = "exp", outcome = c("rmst", "rmtl"), ...)
```

## Arguments

- formula:

  formula with 'Event' outcome

- data:

  data-frame

- model:

  exp ("exp") or identity link ("lin")

- outcome:

  restricted mean time (rmst) or restricted mean time lost (rmtl)

- ...:

  Additional arguments to pass to binregATE

## Details

The first covariate in the specification of the competing risks
regression model must be the treatment effect that is a factor. If the
factor has more than two levels then it uses the mlogit for propensity
score modelling. We consider the outcome mint(T;tau) or
I(epsion==cause1)(t- min(T;t)) that gives years lost due to cause
"cause" depending on the number of causes. The default model is the
exp(X^ beta) and otherwise a linear model is used.

Estimates the ATE using the the standard binary double robust estimating
equations that are IPCW censoring adjusted.

## Author

Thomas Scheike

## Examples

``` r
library(mets); data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
out <- resmeanATE(Event(time,event)~tcell+platelet,data=bmt,time=40,treat.model=tcell~platelet)
summary(out)
#>    n events
#>  408    241
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  2.852563  0.062496  2.730074  2.975052  0.0000
#> tcell1       0.021286  0.122983 -0.219757  0.262329  0.8626
#> platelet     0.303306  0.090772  0.125396  0.481215  0.0008
#> 
#> exp(coeffients):
#>             Estimate     2.5%  97.5%
#> (Intercept) 17.33214 15.33402 19.591
#> tcell1       1.02151  0.80271  1.300
#> platelet     1.35433  1.13360  1.618
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    19.25882  0.95918 17.37887 21.13877  0.0000
#> treat1    19.67316  2.22868 15.30502 24.04129  0.0000
#> treat:1-0  0.41434  2.41151 -4.31213  5.14081  0.8636
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    19.27793  0.95799 17.40030 21.15556  0.0000
#> treat1    20.34004  2.54146 15.35887 25.32122  0.0000
#> treat:1-0  1.06211  2.71020 -4.24979  6.37402  0.6951
#> 
#> 

out1 <- resmeanATE(Event(time,cause)~tcell+platelet,data=bmt,cause=1,time=40,
                   treat.model=tcell~platelet)
summary(out1)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)  2.80626  0.06962  2.66981  2.94271  0.0000
#> tcell1      -0.37413  0.24769 -0.85960  0.11133  0.1309
#> platelet    -0.49164  0.16493 -0.81490 -0.16837  0.0029
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 16.54790 14.43717 18.9672
#> tcell1       0.68788  0.42333  1.1178
#> platelet     0.61162  0.44268  0.8450
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    14.53165  0.95705 12.65587 16.40742  0.0000
#> treat1     9.99609  2.37815  5.33499 14.65718  0.0000
#> treat:1-0 -4.53556  2.57515 -9.58276  0.51164  0.0782
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     14.51355   0.95800  12.63590  16.39120  0.0000
#> treat1      9.36465   2.41708   4.62727  14.10203  0.0001
#> treat:1-0  -5.14890   2.59798 -10.24084  -0.05696  0.0475
#> 
#> 

ratioATE(out,out1,h=function(x) log(x))
#> $ratioG
#>            Estimate Std.Err    2.5%    97.5%   P-value
#> [treat0]    -0.2816  0.1075 -0.4924 -0.07093 9.150e-33
#> [treat1]    -0.6771  0.3183 -1.3008 -0.05328 1.369e-07
#> [treat0].1  -0.3954  0.3363 -1.0545  0.26370 3.333e-05
#> [treat1].1   0.3954  0.3363 -0.2637  1.05454 7.221e-02
#> 
#>  Null Hypothesis: 
#>   [treat0] = 1
#>   [treat1] = 1
#>   [treat0] = 1
#>   [treat1] = 1 
#>  
#> chisq = 170.3366, df = 4, p-value < 2.2e-16
#> 
#> $ratioDR
#>            Estimate Std.Err    2.5%    97.5%   P-value
#> [treat0]    -0.2839  0.1075 -0.4947 -0.07310 7.454e-33
#> [treat1]    -0.7756  0.3460 -1.4538 -0.09753 2.865e-07
#> [treat0].1  -0.4918  0.3619 -1.2011  0.21754 3.756e-05
#> [treat1].1   0.4918  0.3619 -0.2175  1.20109 1.602e-01
#> 
#>  Null Hypothesis: 
#>   [treat0] = 1
#>   [treat1] = 1
#>   [treat0] = 1
#>   [treat1] = 1 
#>  
#> chisq = 168.3765, df = 4, p-value < 2.2e-16
#> 
```
