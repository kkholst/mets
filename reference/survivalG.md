# G-estimator for Cox and Fine-Gray model

Computes G-estimator \$\$ \hat S(t,A=a) = n^{-1} \sum_i \hat
S(t,A=a,Z_i) \$\$ for the Cox model based on phreg og the Fine-Gray
model based on the cifreg function. Gives influence functions of these
risk estimates and SE's are based on these. If first covariate is a
factor then all contrast are computed, and if continuous then considered
covariate values are given by Avalues.

## Usage

``` r
survivalG(
  x,
  data,
  time = NULL,
  Avalues = NULL,
  varname = NULL,
  same.data = TRUE,
  First = FALSE
)
```

## Arguments

- x:

  phreg or cifreg object

- data:

  data frame for risk averaging, must be part of the data used for
  fitting unless same.data=FALSE. When a subset of the data such as the
  treated model should be fitted with cluster(id)

- time:

  for estimate

- Avalues:

  values to compare for first covariate A

- varname:

  if given then averages for this variable, default is first variable

- same.data:

  assumes that same data is used for fitting of survival model and
  averaging.

- First:

  to only use first record for G-averaging, for example when start,stop
  structure is used with phreg

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
bmt$event <- (bmt$cause!=0)*1; bmt$id <- 1:408
dfactor(bmt) <- tcell.f~tcell

fg1 <- cifreg(Event(time,cause)~tcell.f+platelet+age,bmt,cause=1,
              cox.prep=TRUE,propodds=NULL)
summary(survivalG(fg1,bmt,50))
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%   P-value
#> risk0   0.4332 0.02749 0.3793 0.4870 6.316e-56
#> risk1   0.2726 0.05861 0.1577 0.3875 3.297e-06
#> 
#> Average Treatment effect: difference (G-estimator) :
#>     Estimate Std.Err   2.5%    97.5% P-value
#> ps0  -0.1605  0.0635 -0.285 -0.03607 0.01147
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>         Estimate   Std.Err       2.5%       97.5%    P-value
#> [ps0] -0.4630061 0.2211506 -0.8964533 -0.02955888 0.03629353
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 0.6293888 0.4080142 0.9708737 
#> 

ss <- phreg(Surv(time,event)~tcell.f+platelet+age,bmt) 
summary(survivalG(ss,bmt,50))
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%    P-value
#> risk0   0.6539 0.02708 0.6008 0.7069 8.837e-129
#> risk1   0.5639 0.05971 0.4469 0.6810  3.573e-21
#> 
#> Average Treatment effect: difference (G-estimator) :
#>     Estimate Std.Err    2.5%   97.5% P-value
#> ps0 -0.08992 0.06291 -0.2132 0.03338  0.1529
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>         Estimate   Std.Err       2.5%      97.5%   P-value
#> [ps0] -0.1479471 0.1095497 -0.3626606 0.06676628 0.1768548
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 0.8624767 0.6958226 1.0690456 
#> 
#> Average Treatment effect:  survival-difference (G-estimator) :
#>       Estimate   Std.Err        2.5%     97.5%   P-value
#> ps0 0.08992126 0.0629098 -0.03337969 0.2132222 0.1528985
#> 
#> Average Treatment effect: 1-G (survival)-ratio (G-estimator) :
#> log-ratio: 
#>        Estimate   Std.Err        2.5%     97.5%   P-value
#> [ps0] 0.2309406 0.1503721 -0.06378329 0.5256645 0.1245888
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 1.2597844 0.9382083 1.6915826 
#> 

ss <- phreg(Surv(time,event)~strata(tcell.f)+platelet+age,bmt) 
summary(survivalG(ss,bmt,50))
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%    P-value
#> risk0   0.6441 0.02727 0.5906 0.6975 2.397e-123
#> risk1   0.6172 0.07125 0.4776 0.7568  4.611e-18
#> 
#> Average Treatment effect: difference (G-estimator) :
#>     Estimate Std.Err    2.5%  97.5% P-value
#> ps0 -0.02687 0.07622 -0.1763 0.1225  0.7244
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>          Estimate   Std.Err       2.5%     97.5%  P-value
#> [ps0] -0.04261856 0.1228491 -0.2833984 0.1981613 0.728653
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 0.9582769 0.7532197 1.2191590 
#> 
#> Average Treatment effect:  survival-difference (G-estimator) :
#>       Estimate    Std.Err      2.5%    97.5%   P-value
#> ps0 0.02687251 0.07621848 -0.122513 0.176258 0.7244093
#> 
#> Average Treatment effect: 1-G (survival)-ratio (G-estimator) :
#> log-ratio: 
#>         Estimate   Std.Err       2.5%     97.5%   P-value
#> [ps0] 0.07278456 0.2010782 -0.3213214 0.4668906 0.7173734
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 1.0754988 0.7251901 1.5950269 
#> 

sst <- survivalGtime(ss,bmt,n=50)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
plot(sst)


fg1t <- survivalGtime(fg1,bmt,n=50)
plot(fg1t)


#among treated: must specify id to link influence functions
ss <- phreg(Surv(time,event)~tcell.f+platelet+age+cluster(id),bmt) 
summary(survivalG(ss,subset(bmt,tcell==1),50))
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%   P-value
#> risk0   0.6662 0.03407 0.5995 0.7330 3.661e-85
#> risk1   0.5749 0.05749 0.4622 0.6876 1.518e-23
#> 
#> Average Treatment effect: difference (G-estimator) :
#>     Estimate Std.Err    2.5%   97.5% P-value
#> ps0 -0.09134 0.06417 -0.2171 0.03442  0.1546
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>         Estimate   Std.Err       2.5%      97.5%   P-value
#> [ps0] -0.1474634 0.1081918 -0.3595155 0.06458876 0.1728887
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 0.8628941 0.6980145 1.0667203 
#> 
#> Average Treatment effect:  survival-difference (G-estimator) :
#>       Estimate    Std.Err        2.5%     97.5%   P-value
#> ps0 0.09134406 0.06416628 -0.03441954 0.2171077 0.1545761
#> 
#> Average Treatment effect: 1-G (survival)-ratio (G-estimator) :
#> log-ratio: 
#>        Estimate   Std.Err        2.5%     97.5%   P-value
#> [ps0] 0.2419052 0.1620237 -0.07565548 0.5594659 0.1354311
#> ratio: 
#>  Estimate      2.5%     97.5% 
#> 1.2736735 0.9271356 1.7497378 
#> 
```
