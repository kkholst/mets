# Analysis of multivariate binomial data: family analysis

## Overview

When looking at multivariate binomial data with the aim of learning
about the dependence that is present, possibly after correcting for some
covariates many models and methods are available:

- Random-effects models logistic regression covered elsewhere (glmer in
  lme4).

in the mets package you can fit the

- Pairwise odds ratio model

- Bivariate Probit model

  - With random effects
  - Special functionality for polygenic random effects modelling such as
    ACE, ADE ,AE and so forth.

- Additive gamma random effects model

  - Special functionality for polygenic random effects modelling such as
    ACE, ADE ,AE and so forth.

These last three models are all fitted in the mets package using
composite likelihoods based on pairs within clusters. The models can be
fitted specifically based on specifying which pairs one wants to use for
the composite score.

The models are described in futher details in the binomial-twin
vignette.

## Simulated family data

We start by simulating family data with and additive gamma structure on
ACE form. Here 10000 families consisting of two parents and two
children. The response is ybin and there is one covariate x.

``` r
 library(mets)
 library(timereg)
#> Loading required package: survival
#> 
#> Attaching package: 'timereg'
#> The following objects are masked from 'package:mets':
#> 
#>     Event, event.split, kmplot, plotConfregion
 set.seed(100)
 data <- simbinClaytonOakes.family.ace(500,2,1,beta=NULL,alpha=NULL)
 data$number <- c(1,2,3,4)
 data$child <- 1*(data$number==3)
 head(data)
#>   ybin x   type cluster number child
#> 1    1 1 mother       1      1     0
#> 2    1 0 father       1      2     0
#> 3    1 0  child       1      3     1
#> 4    1 1  child       1      4     0
#> 5    1 1 mother       2      1     0
#> 6    1 0 father       2      2     0
```

We fit the marginal models, and here find a covariate effect at 0.3 for
x. The marginals can be specified as desired.

``` r
 aa <- margbin <- glm(ybin~x,data=data,family=binomial())
 summary(aa)
#> 
#> Call:
#> glm(formula = ybin ~ x, family = binomial(), data = data)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.49908    0.06387   7.814 5.52e-15 ***
#> x            0.17531    0.09355   1.874   0.0609 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 2610.2  on 1999  degrees of freedom
#> Residual deviance: 2606.7  on 1998  degrees of freedom
#> AIC: 2610.7
#> 
#> Number of Fisher Scoring iterations: 4
```

## Additive Gamma model

For the Additive Gamma model we set-up the random effects included in
such a family to make the ACE valid using some special functions for
this.

The model is constructed with one enviromental effect shared by all in
the family and 8 genetic random effects with size 1/4 genetic variance.
Looking at the first family we see that the mother and father both share
half the genes with the children and that the two children also share
half their genes with this specification. Below we also show an
alternative specification of this model using all pairs.

``` r
# make ace random effects design
out <- ace.family.design(data,member="type",id="cluster")
out$pardes
#>       [,1] [,2]
#>  [1,] 0.25    0
#>  [2,] 0.25    0
#>  [3,] 0.25    0
#>  [4,] 0.25    0
#>  [5,] 0.25    0
#>  [6,] 0.25    0
#>  [7,] 0.25    0
#>  [8,] 0.25    0
#>  [9,] 0.00    1
head(out$des.rv,4)
#>      m1 m2 m3 m4 f1 f2 f3 f4 env
#> [1,]  1  1  1  1  0  0  0  0   1
#> [2,]  0  0  0  0  1  1  1  1   1
#> [3,]  1  1  0  0  1  1  0  0   1
#> [4,]  1  0  1  0  1  0  1  0   1
```

We can now fit the model calling the two-stage function

``` r
# fitting ace model for family structure
ts <- binomial.twostage(margbin,data=data,clusters=data$cluster,
theta=c(2,1),random.design=out$des.rv,theta.des=out$pardes)
summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                theta        se
#> dependence1 2.817486 0.8981749
#> dependence2 1.137637 0.2806719
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err   2.5%  97.5%   P-value
#> dependence1   0.7124 0.09262 0.5308 0.8939 1.457e-14
#> dependence2   0.2876 0.09262 0.1061 0.4692 1.899e-03
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    3.955  0.8668 2.256 5.654 5.051e-06
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
# true variance parameters
c(2,1)
#> [1] 2 1
# total variance 
3
#> [1] 3
```

## Pairwise fitting

We now specify the same model via extracting all pairs.  
The random effecs structure is typically simpler when just looking at
pairs, but we start by specifying the random effects jointly for whole
family. A special function writes up all combinations of pairs. There
are 6 pairs within each family, and we keep track of who belongs to the
different families. We first simply give the pairs and we then should
get the same result as before.

``` r
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
#>  [1] 1 2 1 3 1 4 2 3 2 4 3 4 5 6 5 7 5 8 6 7
pairs <- mm$pairs
dim(pairs)
#> [1] 3000    2
head(pairs,12)
#>       [,1] [,2]
#>  [1,]    1    2
#>  [2,]    1    3
#>  [3,]    1    4
#>  [4,]    2    3
#>  [5,]    2    4
#>  [6,]    3    4
#>  [7,]    5    6
#>  [8,]    5    7
#>  [9,]    5    8
#> [10,]    6    7
#> [11,]    6    8
#> [12,]    7    8
```

Now with the pairs we fit the model

``` r
tsp <- binomial.twostage(margbin,data=data,clusters=data$cluster,theta=c(2,1),detail=0,
        random.design=out$des.rv,theta.des=out$pardes,pairs=pairs)
summary(tsp)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                theta        se
#> dependence1 2.817486 0.8981749
#> dependence2 1.137637 0.2806719
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err   2.5%  97.5%   P-value
#> dependence1   0.7124 0.09262 0.5308 0.8939 1.457e-14
#> dependence2   0.2876 0.09262 0.1061 0.4692 1.899e-03
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    3.955  0.8668 2.256 5.654 5.051e-06
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Here a random sample of pairs are given instead and we get other
estimates.

``` r
set.seed(100)
ssid <- sort(sample(1:nrow(pairs),nrow(pairs)/2))
tsd <- binomial.twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1),step=1.0,
               random.design=out$des.rv,iid=1,Nit=10,
               theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 theta        se
#> dependence1 3.0275402 1.6045129
#> dependence2 0.8135103 0.3895686
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5%   P-value
#> dependence1   0.7882  0.1511  0.49210 1.0843 1.816e-07
#> dependence2   0.2118  0.1511 -0.08431 0.5079 1.609e-01
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5% P-value
#> p1    3.841   1.402 1.093 6.589 0.00615
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

To specify such a model when only the pairs are available we show how to
specify the model. We here use the same marginal “aa” to make the
results comparable. The marginal can also be fitted based on available
data.

We start by selecting the data related to the pairs, and sets up new
id’s and to start we specify the model using the full design with 9
random effects. Below we show how one can use with only the random
effects needed for each pair, which is typically simpler.

``` r
head(pairs[ssid,])
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    2    3
#> [4,]    2    4
#> [5,]    3    4
#> [6,]    6    8
ids <- sort(unique(c(pairs[ssid,])))

pairsids <- c(pairs[ssid,])
pair.new <- matrix(fast.approx(ids,c(pairs[ssid,])),ncol=2)
head(pair.new)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    2    3
#> [4,]    2    4
#> [5,]    3    4
#> [6,]    5    7

dataid <- dsort(data[ids,],"cluster")
outid <- ace.family.design(dataid,member="type",id="cluster")
outid$pardes
#>       [,1] [,2]
#>  [1,] 0.25    0
#>  [2,] 0.25    0
#>  [3,] 0.25    0
#>  [4,] 0.25    0
#>  [5,] 0.25    0
#>  [6,] 0.25    0
#>  [7,] 0.25    0
#>  [8,] 0.25    0
#>  [9,] 0.00    1
head(outid$des.rv)
#>      m1 m2 m3 m4 f1 f2 f3 f4 env
#> [1,]  1  1  1  1  0  0  0  0   1
#> [2,]  0  0  0  0  1  1  1  1   1
#> [3,]  1  1  0  0  1  1  0  0   1
#> [4,]  1  0  1  0  1  0  1  0   1
#> [5,]  0  0  0  0  1  1  1  1   1
#> [6,]  1  1  0  0  1  1  0  0   1
```

Now fitting the model with the data set up

``` r
aa <- glm(ybin~x,data=dataid,family=binomial())
tsdid <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
         theta=c(2,1),random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 theta       se
#> dependence1 3.1389859 1.649758
#> dependence2 0.8010544 0.399870
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5%   P-value
#> dependence1   0.7967  0.1492  0.50421 1.0892 9.358e-08
#> dependence2   0.2033  0.1492 -0.08917 0.4958 1.731e-01
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%  P-value
#> p1     3.94   1.438 1.121 6.759 0.006153
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

We now specify the design specifically using the pairs. The
random.design and design on the parameters are now given for each pair,
as a 3 dimensional matrix. with a direct specification of random.design
and the design on the parameters theta.design. In addition we need also
to give the number of random effects for each pair. These basic things
are constructed by certain functions for the ACE design.

``` r
pair.types <-  matrix(dataid[c(t(pair.new)),"type"],byrow=T,ncol=2)
head(pair.new,7)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    2    3
#> [4,]    2    4
#> [5,]    3    4
#> [6,]    5    7
#> [7,]    6    7
head(pair.types,7)
#>      [,1]     [,2]    
#> [1,] "mother" "father"
#> [2,] "mother" "child" 
#> [3,] "father" "child" 
#> [4,] "father" "child" 
#> [5,] "child"  "child" 
#> [6,] "father" "child" 
#> [7,] "child"  "child"
###
theta.des  <- rbind( c(rbind(c(1,0),  c(1,0),  c(0,1),  c(0,0))),
             c(rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))))
random.des <- rbind( 
        c(1,0,1,0),c(0,1,1,0),
        c(1,1,0,1),c(1,0,1,1))
mf <- 1*(pair.types[,1]=="mother" & pair.types[,2]=="father")
##          pair, rv related to pairs,  theta.des related to pair 
pairs.new <- cbind(pair.new,(mf==1)*1+(mf==0)*3,(mf==1)*2+(mf==0)*4,(mf==1)*1+(mf==0)*2,(mf==1)*3+(mf==0)*4)
```

pairs.new is matrix with

- columns 1:2 giving the indeces of the data points

- columns 3:4 giving the indeces of the random.design for the different
  pairs

- columns 5 giving the indeces of the theta.des written as rows

- columns 6 giving the number of random variables for this pair

The length of all rows of theta.des are the maximum number of random
effects $\times$ the number of parameters. These two numbers are given
in the call. In this case 4 $\times$ 2. So theta.des has rows of length
$8$, possibly including some 0’s for rows not relevant due to fewer
random effects, as is the case here for pairs that do not share genetic
effects.

For pair 1 that is a mother/farther pair, we see that they share 1
environmental random effect of size 1. There are also two genetic
effects that are unshared between the two. So a total of 3 random
effects are needed here. The theta.des relates the 3 random effects to
possible relationships in the parameters. Here the genetic effects are
full and so is the environmental effect. In contrast we also consider a
mother/child pair that share half the genes, now with random effects
with (1/2) gene variance. We there need 4 random effects, 2 non-shared
half-gene, 1 shared half-gene, and one shared full environmental effect.

``` r
# 3 rvs here 
random.des
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    1    0
#> [2,]    0    1    1    0
#> [3,]    1    1    0    1
#> [4,]    1    0    1    1
theta.des
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]  1.0  1.0  0.0    0    0    0    1    0
#> [2,]  0.5  0.5  0.5    0    0    0    0    1

head(pairs.new)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    2    1    2    1    3
#> [2,]    1    3    3    4    2    4
#> [3,]    2    3    3    4    2    4
#> [4,]    2    4    3    4    2    4
#> [5,]    3    4    3    4    2    4
#> [6,]    5    7    3    4    2    4
```

Now fitting the model, and we see that it is a lot quicker due to the
fewer random effects needed for pairs. We need to also specify the
number of parameters in this case.

``` r
tsdid2 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster, theta=c(2,1),
           random.design=random.des,theta.des=theta.des,pairs=pairs.new,dim.theta=2)
summary(tsdid2)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 theta       se
#> dependence1 3.1389859 1.649758
#> dependence2 0.8010544 0.399870
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5%   P-value
#> dependence1   0.7967  0.1492  0.50421 1.0892 9.358e-08
#> dependence2   0.2033  0.1492 -0.08917 0.4958 1.731e-01
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%  P-value
#> p1     3.94   1.438 1.121 6.759 0.006153
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

The same model can be specifed even simpler via the kinship coefficient.
For this speicification there are 4 random effects for each pair, but
some have variance 0. The mother-father pair, here shares a random
effect with variance 0, and have two non-shared genetic effects with
full variance, in addition to a fully shared environmental effect.

``` r
kinship  <- rep(0.5,nrow(pair.types))
kinship[pair.types[,1]=="mother" & pair.types[,2]=="father"] <- 0
head(kinship,n=10)
#>  [1] 0.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5

out <- make.pairwise.design(pair.new,kinship,type="ace") 
```

``` r
tsdid3 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
             theta=c(2,1)/9,random.design=out$random.design,
             theta.des=out$theta.des,pairs=out$new.pairs,dim.theta=2)
summary(tsdid3)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 theta       se
#> dependence1 3.1389859 1.649758
#> dependence2 0.8010544 0.399870
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5%   P-value
#> dependence1   0.7967  0.1492  0.50421 1.0892 9.358e-08
#> dependence2   0.2033  0.1492 -0.08917 0.4958 1.731e-01
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%  P-value
#> p1     3.94   1.438 1.121 6.759 0.006153
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

### Pairwise dependence modelling

``` r
library(mets)
set.seed(1000)
data <- simbinClaytonOakes.family.ace(500,2,1,beta=NULL,alpha=NULL)
head(data)
#>   ybin x   type cluster
#> 1    0 1 mother       1
#> 2    0 0 father       1
#> 3    0 1  child       1
#> 4    0 0  child       1
#> 5    1 0 mother       2
#> 6    1 0 father       2
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
#>  [1] 1 2 1 3 1 4 2 3 2 4 3 4 5 6 5 7 5 8 6 7
pairs <- mm$pairs
dim(pairs)
#> [1] 3000    2
head(pairs,12)
#>       [,1] [,2]
#>  [1,]    1    2
#>  [2,]    1    3
#>  [3,]    1    4
#>  [4,]    2    3
#>  [5,]    2    4
#>  [6,]    3    4
#>  [7,]    5    6
#>  [8,]    5    7
#>  [9,]    5    8
#> [10,]    6    7
#> [11,]    6    8
#> [12,]    7    8
```

``` r
 dtypes <- interaction( data[pairs[,1],"type"], data[pairs[,2],"type"])
 dtypes <- droplevels(dtypes)
 table(dtypes)
#> dtypes
#>   child.child  father.child  mother.child mother.father 
#>           500          1000          1000           500
 dm <- model.matrix(~-1+factor(dtypes))
```

Now with the pairs we fit the model

``` r
aa <- glm(ybin~x,data=data,family=binomial())

tsp <- binomial.twostage(aa,data=data, clusters=data$cluster,
         theta.des=dm,pairs=cbind(pairs,1:nrow(dm)))
summary(tsp)
#> Dependence parameter for Odds-Ratio (Plackett) model
#> $estimates
#>                                theta        se
#> factor(dtypes)child.child   3.520109 0.7143594
#> factor(dtypes)father.child  3.871698 0.6055928
#> factor(dtypes)mother.child  3.922180 0.6233304
#> factor(dtypes)mother.father 1.328730 0.2657202
#> 
#> $log.or
#>                             Estimate Std.Err    2.5%  97.5%   P-value
#> factor(dtypes)child.child     1.2585  0.2029  0.8607 1.6562 5.596e-10
#> factor(dtypes)father.....     1.3537  0.1564  1.0471 1.6603 4.952e-18
#> factor(dtypes)mother.....     1.3666  0.1589  1.0552 1.6781 8.017e-18
#> factor(dtypes)mother......1   0.2842  0.2000 -0.1077 0.6762 1.552e-01
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

## Pairwise odds ratio model

To fit the pairwise odds-ratio model in the case of a pair-specification
there are two options for fitting the model.

1.  One option is to set up some artificial data similar to twin data
    with

- a pair-cluster-id (clusters)
- with a cluster-id to get GEE type standard errors (se.cluster)

2.  We can also use the specify the design via the theta.des that is
    also a matrix of dimension pairs x design with the design for POR
    model.

Starting by the second option. We need to start by specify the design of
the odds-ratio of each pair. We set up the data and find all
combinations within the pairs. Subsequently, we remove all the empty
groups, by grouping together the factor levels 4:9, and then we
construct the design.

``` r
tdp <-cbind( dataid[pair.new[,1],],dataid[pair.new[,2],])
names(tdp) <- c(paste(names(dataid),"1",sep=""),
        paste(names(dataid),"2",sep=""))
tdp <-transform(tdp,tt=interaction(type1,type2))
dlevel(tdp)
#> tt #levels=:6 
#> [1] "child.child"   "father.child"  "mother.child"  "child.father" 
#> [5] "father.father" "mother.father"
#> -----------------------------------------
drelevel(tdp,newlevels=list(mother.father=4:9)) <-  obs.types~tt
dtable(tdp,~tt+obs.types)
#> 
#>               obs.types mother.father child.child father.child mother.child
#> tt                                                                         
#> child.child                         0         232            0            0
#> father.child                        0           0          525            0
#> mother.child                        0           0            0          510
#> child.father                        0           0            0            0
#> father.father                       0           0            0            0
#> mother.father                     233           0            0            0
tdp <- model.matrix(~-1+factor(obs.types),tdp)
```

We then can fit the pairwise model using the pairs and the pair-design
for descrbing the OR. The results are consistent with the the ACE model
as the mother-father have a lower dependence as is due only the
environmental effects. All other combinations should have the same
dependence as also seem to be the case.

To fit the OR model it is generally recommended to use the var.link to
use the parmetrization with log-odd-ratio regression.

``` r
###porpair <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
###           theta.des=tdp,pairs=pair.new,model="or",var.link=1)
###summary(porpair)
```

## SessionInfo

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] timereg_2.0.7  survival_3.8-3 mets_1.3.9    
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.5           knitr_1.51          rlang_1.1.6        
#>  [4] xfun_0.55           textshaping_1.0.4   jsonlite_2.0.0     
#>  [7] listenv_0.10.0      future.apply_1.20.1 lava_1.8.2         
#> [10] htmltools_0.5.9     ragg_1.5.0          sass_0.4.10        
#> [13] rmarkdown_2.30      grid_4.5.2          evaluate_1.0.5     
#> [16] jquerylib_0.1.4     fastmap_1.2.0       numDeriv_2016.8-1.1
#> [19] yaml_2.3.12         mvtnorm_1.3-3       lifecycle_1.0.4    
#> [22] compiler_4.5.2      codetools_0.2-20    fs_1.6.6           
#> [25] htmlwidgets_1.6.4   Rcpp_1.1.0          future_1.68.0      
#> [28] lattice_0.22-7      systemfonts_1.3.1   digest_0.6.39      
#> [31] R6_2.6.1            parallelly_1.46.0   parallel_4.5.2     
#> [34] splines_4.5.2       Matrix_1.7-4        bslib_0.9.0        
#> [37] tools_4.5.2         globals_0.18.0      pkgdown_2.2.0      
#> [40] cachem_1.1.0        desc_1.4.3
```
