# Fits two-stage binomial for describing depdendence in binomial data using marginals that are on logistic form using the binomial.twostage funcion, but call is different and easier and the data manipulation is build into the function. Useful in particular for family design data.

If clusters contain more than two times, the algoritm uses a compososite
likelihood based on the pairwise bivariate models.

## Usage

``` r
easy.binomial.twostage(
  margbin = NULL,
  data = parent.frame(),
  method = "nr",
  response = "response",
  id = "id",
  Nit = 60,
  detail = 0,
  silent = 1,
  weights = NULL,
  control = list(),
  theta = NULL,
  theta.formula = NULL,
  desnames = NULL,
  deshelp = 0,
  var.link = 1,
  iid = 1,
  step = 1,
  model = "plackett",
  marginal.p = NULL,
  strata = NULL,
  max.clust = NULL,
  se.clusters = NULL
)
```

## Arguments

- margbin:

  Marginal binomial model

- data:

  data frame

- method:

  Scoring method

- response:

  name of response variable in data frame

- id:

  name of cluster variable in data frame

- Nit:

  Number of iterations

- detail:

  Detail for more output for iterations

- silent:

  Debug information

- weights:

  Weights for log-likelihood, can be used for each type of outcome in
  2x2 tables.

- control:

  Optimization arguments

- theta:

  Starting values for variance components

- theta.formula:

  design for depedence, either formula or design function

- desnames:

  names for dependence parameters

- deshelp:

  if 1 then prints out some data sets that are used, on on which the
  design function operates

- var.link:

  Link function for variance

- iid:

  Calculate i.i.d. decomposition

- step:

  Step size

- model:

  model

- marginal.p:

  vector of marginal probabilities

- strata:

  strata for fitting

- max.clust:

  max clusters used for i.i.d. decompostion

- se.clusters:

  clusters for iid decomposition for roubst standard errors

## Details

The reported standard errors are based on the estimated information from
the likelihood assuming that the marginals are known. This gives correct
standard errors in the case of the plackett distribution (OR model for
dependence), but incorrect for the clayton-oakes types model. The OR
model is often known as the ALR model. Our fitting procedures gives
correct standard errors due to the ortogonality and is fast.

## Examples

``` r
data(twinstut)
twinstut0 <- subset(twinstut, tvparnr<4000)
twinstut <- twinstut0
twinstut$binstut <- (twinstut$stutter=="yes")*1
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
bin <- binomial.twostage(margbin,data=twinstut,var.link=1,
             clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
                   method="nr")
summary(bin)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                    theta        se
#> factor(zyg)dz -0.2853738 0.9894082
#> factor(zyg)mz  3.3391390 0.5590195
#> factor(zyg)os  0.4920396 0.7634939
#> 
#> $or
#>               Estimate Std.Err   2.5%  97.5% P-value
#> factor(zyg)dz   0.7517  0.7438 -0.706  2.209 0.31216
#> factor(zyg)mz  28.1948 15.7615 -2.697 59.087 0.07364
#> factor(zyg)os   1.6356  1.2488 -0.812  4.083 0.19027
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
lava::estimate(coef=bin$theta,vcov=bin$var.theta,f=function(p) exp(p))
#>      Estimate Std.Err   2.5%  97.5% P-value
#> NA     0.7517  0.7438 -0.706  2.209 0.31216
#> NA.1  28.1948 15.7615 -2.697 59.087 0.07364
#> NA.2   1.6356  1.2488 -0.812  4.083 0.19027

twinstut$cage <- scale(twinstut$age)
theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
             clusters=twinstut$tvparnr,theta.des=theta.des,detail=0)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                    theta        se
#> factor(zyg)dz -0.2684851 0.9930894
#> factor(zyg)mz  3.4239727 0.5773886
#> factor(zyg)os  0.4778091 0.7628390
#> cage           0.2519096 0.4821619
#> 
#> $or
#>               Estimate Std.Err     2.5%  97.5% P-value
#> factor(zyg)dz   0.7645  0.7593 -0.72357  2.253 0.31395
#> factor(zyg)mz  30.6911 17.7207 -4.04082 65.423 0.08328
#> factor(zyg)os   1.6125  1.2301 -0.79843  4.024 0.18989
#> cage            1.2865  0.6203  0.07073  2.502 0.03808
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
             clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                          theta        se
#> factor(zyg)dz      -0.27701246 1.0713974
#> factor(zyg)mz       3.52484148 0.5951743
#> factor(zyg)os       0.48859941 0.7644050
#> cage                0.06420907 3.6225641
#> factor(zyg)mz:cage  0.49441325 3.6865646
#> factor(zyg)os:cage -0.12312666 3.6921748
#> 
#> $or
#>                    Estimate Std.Err     2.5%  97.5% P-value
#> factor(zyg)dz        0.7580  0.8122  -0.8338  2.350 0.35063
#> factor(zyg)mz       33.9484 20.2052  -5.6531 73.550 0.09292
#> factor(zyg)os        1.6300  1.2460  -0.8121  4.072 0.19080
#> cage                 1.0663  3.8628  -6.5046  8.637 0.78251
#> factor(zyg)mz:cage   1.6395  6.0443 -10.2070 13.486 0.78619
#> factor(zyg)os:cage   0.8842  3.2644  -5.5140  7.282 0.78651
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,
                              response="binstut",id="tvparnr",var.link=1,
                theta.formula=~-1+factor(zyg1))
summary(out)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                     theta        se
#> factor(zyg1)dz -0.2853738 0.9894082
#> factor(zyg1)mz  3.3391390 0.5590195
#> factor(zyg1)os  0.4920396 0.7634939
#> 
#> $or
#>                Estimate Std.Err   2.5%  97.5% P-value
#> factor(zyg1)dz   0.7517  0.7438 -0.706  2.209 0.31216
#> factor(zyg1)mz  28.1948 15.7615 -2.697 59.087 0.07364
#> factor(zyg1)os   1.6356  1.2488 -0.812  4.083 0.19027
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

## refers to zygosity of first subject in eash pair : zyg1
## could also use zyg2 (since zyg2=zyg1 within twinpair's))
## do not run t save time
# desfs <- function(x,num1="zyg1",namesdes=c("mz","dz","os"))
#     c(x[num1]=="mz",x[num1]=="dz",x[num1]=="os")*1
#
#out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
#                               data=twinstut, response="binstut",id="tvparnr",
#                               var.link=1,theta.formula=desfs,
#                               desnames=c("mz","dz","os"))
#summary(out3)

 ## Reduce Ex.Timings
n <- 1000
set.seed(100)
dd <- simBinFam(n,beta=0.3)
binfam <- fast.reshape(dd,varying=c("age","x","y"))
## mother, father, children  (ordered)
head(binfam)
#>   id      age x y num
#> 1  1 26.68147 1 1   m
#> 2  1 32.93535 0 1   f
#> 3  1  9.09635 0 1  b1
#> 4  1 11.05365 1 1  b2
#> 5  2 26.06935 1 1   m
#> 6  2 28.33452 0 1   f

########### ########### ########### ########### ########### ###########
####  simple analyses of binomial family data
########### ########### ########### ########### ########### ###########
desfs <- function(x,num1="num1",num2="num2")
{
     pp <- 1*(((x[num1]=="m")*(x[num2]=="f"))|(x[num1]=="f")*(x[num2]=="m"))
     pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
     cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
     c(pp,pc,cc)
}

ud <- easy.binomial.twostage(y~+1,data=binfam,
     response="y",id="id",
     theta.formula=desfs,desnames=c("pp","pc","cc"))
summary(ud)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>        theta        se
#> pp 0.6561647 0.2508517
#> pc 0.5022035 0.1393882
#> cc 0.8892008 0.2124875
#> 
#> $or
#>    Estimate Std.Err   2.5% 97.5%   P-value
#> pp    1.927  0.4835 0.9798 2.875 6.708e-05
#> pc    1.652  0.2303 1.2009 2.104 7.273e-13
#> cc    2.433  0.5170 1.4198 3.447 2.524e-06
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

udx <- easy.binomial.twostage(y~+x,data=binfam,
     response="y",id="id",
     theta.formula=desfs,desnames=c("pp","pc","cc"))
summary(udx)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>        theta        se
#> pp 0.5704932 0.2541853
#> pc 0.5125353 0.1410238
#> cc 0.8840448 0.2153984
#> 
#> $or
#>    Estimate Std.Err   2.5% 97.5%   P-value
#> pp    1.769  0.4497 0.8878 2.651 8.350e-05
#> pc    1.670  0.2354 1.2081 2.131 1.331e-12
#> cc    2.421  0.5214 1.3987 3.443 3.441e-06
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

########### ########### ########### ########### ########### ###########
####  now allowing parent child POR to be different for mother and father
########### ########### ########### ########### ########### ###########

desfsi <- function(x,num1="num1",num2="num2")
{
    pp <- (x[num1]=="m")*(x[num2]=="f")*1
    mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
    fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
    cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
    c(pp,mc,fc,cc)
}

udi <- easy.binomial.twostage(y~+1,data=binfam,
     response="y",id="id",
     theta.formula=desfsi,desnames=c("pp","mother-child","father-child","cc"))
summary(udi)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                  theta        se
#> pp           0.6561647 0.2508517
#> mother-child 0.7897223 0.1715685
#> father-child 0.1784171 0.1944482
#> cc           0.8892008 0.2124875
#> 
#> $or
#>              Estimate Std.Err   2.5% 97.5%   P-value
#> pp              1.927  0.4835 0.9798 2.875 6.708e-05
#> mother-child    2.203  0.3779 1.4621 2.944 5.590e-09
#> father-child    1.195  0.2324 0.7398 1.651 2.707e-07
#> cc              2.433  0.5170 1.4198 3.447 2.524e-06
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

##now looking to see if interactions with age or age influences marginal models
##converting factors to numeric to make all involved covariates numeric
##to use desfai2 rather then desfai that works on binfam

nbinfam <- binfam
nbinfam$num <- as.numeric(binfam$num)
head(nbinfam)
#>   id      age x y num
#> 1  1 26.68147 1 1   1
#> 2  1 32.93535 0 1   2
#> 3  1  9.09635 0 1   3
#> 4  1 11.05365 1 1   4
#> 5  2 26.06935 1 1   1
#> 6  2 28.33452 0 1   2

desfsai <- function(x,num1="num1",num2="num2")
{
    pp <- (x[num1]=="m")*(x[num2]=="f")*1
### av age for pp=1 i.e parent pairs
    agepp <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-30)*pp
    mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
    fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
    cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
    agecc <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-12)*cc
    c(pp,agepp,mc,fc,cc,agecc)
}

desfsai2 <- function(x,num1="num1",num2="num2")
{
    pp <- (x[num1]==1)*(x[num2]==2)*1
    agepp <- (((x["age1"]+x["age2"]))/2-30)*pp ### av age for pp=1 i.e parent pairs
    mc <- (x[num1]==1)*(x[num2]==3 | x[num2]==4)*1
    fc <- (x[num1]==2)*(x[num2]==3 | x[num2]==4)*1
    cc <- (x[num1]==3)*(x[num2]==3 | x[num2]==4)*1
    agecc <- ((x["age1"]+x["age2"])/2-12)*cc ### av age for children
    c(pp,agepp,mc,fc,cc,agecc)
}

udxai2 <- easy.binomial.twostage(y~+x+age,data=binfam,
     response="y",id="id",
     theta.formula=desfsai,
     desnames=c("pp","pp-age","mother-child","father-child","cc","cc-age"))
summary(udxai2)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                    theta         se
#> pp            0.38859947 0.40359602
#> pp-age       -0.05025909 0.08249533
#> mother-child  0.80559724 0.17583430
#> father-child  0.18500933 0.19599386
#> cc            0.85277851 0.22046314
#> cc-age       -0.04548620 0.07212976
#> 
#> $or
#>              Estimate Std.Err   2.5% 97.5%   P-value
#> pp             1.4749 0.59527 0.3082 2.642 1.322e-02
#> pp-age         0.9510 0.07845 0.7972 1.105 8.086e-34
#> mother-child   2.2380 0.39352 1.4667 3.009 1.292e-08
#> father-child   1.2032 0.23583 0.7410 1.665 3.357e-07
#> cc             2.3462 0.51724 1.3324 3.360 5.736e-06
#> cc-age         0.9555 0.06892 0.8204 1.091 1.048e-43
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```
