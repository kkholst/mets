# Analysis of bivariate binomial data: Twin analysis

## Overview

When looking at bivariate binomial data with the aim of learning about
the dependence that is present, possibly after correcting for some
covariates many models are available.

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

Typically it can be hard or impossible to specify random effects models
with special structure among the parameters of the random effects. This
is possible in our models.

To be concrete about the model structure assume that we have paired
binomial data $(Y_{1},Y_{2},X_{1},X_{)}2$ where the responses are
$Y_{1},Y_{2}$ and we have covariates $X_{1},X_{2}$.

We start by giving a brief description of these different models. First
we for bivariate data one can specify the marginal probability using
logistic regression models
$$logit\left( P\left( Y_{i} = 1|X_{i} \right) \right) = \alpha_{i} + X_{i}^{T}\beta i = 1,2.$$
These model can be estimated under working independence .

A typical twin analysis will typically consist of looking at both

- Pairwise odds ratio model

- Bivariate Probit model

The additive gamma can be used for the same as the bivariate probit
model but is more restrictive in terms of dependence structure, but is
nevertheless still valuable to have also as a check of results of the
bivariate probit model.

## Biprobit with random effects

For these model we assume that given random effects $Z$ and a covariate
vector $V_{12}$ we have independent logistic regression models
$$probit\left( P\left( Y_{i} = 1|X_{i},Z \right) \right) = \alpha_{i} + X_{i}^{T}\beta + V_{12}^{T}Zi = 1,2.$$
where Z\$ is a bivariate normal distribution with some covariance
$\Sigma$. The general covariance structure $\Sigma$ makes the model very
flexible.

We note that

- Paramters $\beta$ are subject specific
- The $\Sigma$ will reflect dependence

The more standard link function $logit$ rather than the $probit$ link is
often used and implemented in for example . The advantage is that one
now gets an odds-ratio interpretation of the subject specific effects,
but one then needs numerical integration to fit the model.

### Pairwise odds ratio model

Now the pairwise odds ratio model the specifies that given \$ X_1, X_2
\$ the marginal models are
$$logit\left( P\left( Y_{i} = 1|X_{i} \right) \right) = \alpha_{i} + X_{i}^{T}\beta i = 1,2$$

The primary object of interest are the odds ratio between $Y_{1}$ and
$Y_{2}$$$\gamma_{12} = \frac{P\left( Y_{ki} = 1,Y_{kj} = 1 \right)P\left( Y_{ki} = 0,Y_{kj} = 0 \right)}{P\left( Y_{ki} = 1,Y_{kj} = 0 \right)P\left( Y_{ki} = 0,Y_{kj} = 1 \right)}$$
given $X_{ki}$, $X_{kj}$, and $Z_{kji}$.

We model the odds ratio with the regression
$$\gamma_{12} = \exp\left( Z_{12}^{T}\lambda \right)$$ Where $Z_{12}$
are some covarites that may influence the odds-ratio between between
$Y_{1}$ and $Y_{2}$ and contains the marginal covariates, . This
odds-ratio is given covariates as well as marginal covariates. The
odds-ratio and marginals specify the joint bivariate distribution via
the so-called Placckett-distribution.

One way of fitting this model is the ALR algoritm, the alternating
logistic regression ahd this has been described in several papers . We
here simply estimate the parameters in a two stage-procedure

- Estimating the marginal parameters via GEE
- Using marginal estimates, estimate dependence parameters

This gives efficient estimates of the dependence parameters because of
orthogonality, but some efficiency may be gained for the marginal
parameters by using the full likelihood or iterative fitting such as for
the ALR.

The pairwise odds-ratio model is very useful, but one do not have a
random effects model.

### Additive gamma model

Again we operate under marginal logistic regression models are
$$logit\left( P\left( Y_{i} = 1|X_{i} \right) \right) = \alpha_{i} + X_{i}^{T}\beta i = 1,2$$

First with just one random effect $Z$ we assume that conditional on $Z$
the responses are independent and follow the model
$$logit\left( P\left( Y_{i} = 1|X_{i},Z \right) \right) = exp\left( - Z \cdot \Psi^{- 1}\left( \lambda_{\bullet},\lambda_{\bullet},P\left( Y_{i} = 1|X_{i} \right) \right) \right)$$
where $\Psi$ is the laplace transform of $Z$ where we assume that $Z$ is
gamma distributed with variance $\lambda_{\bullet}^{- 1}$ and mean 1. In
general $\Psi\left( \lambda_{1},\lambda_{2} \right)$ is the laplace
transform of a Gamma distributed random effect with $Z$ with mean
$\lambda_{1}/\lambda_{2}$ and variance $\lambda_{1}/\lambda_{2}^{2}$.

We fit this model by

- Estimating the marginal parameters via GEE
- Using marginal estimates, estimate dependence parameters

To deal with multiple random effects we consider random effects
$Z_{i}i = 1,...,d$ such that $Z_{i}$ is gamma distributed with \* mean:
$\lambda_{j}/\lambda_{\bullet}$ \* variance:
$\lambda_{j}/\lambda_{\bullet}^{2}$, where we define the scalar
$\lambda_{\bullet}$ below.

Now given a cluster-specific design vector $V_{12}$ we assume that
$$V_{12}^{T}Z$$ is gamma distributed with mean 1 and variance
$\lambda_{\bullet}^{- 1}$ such that critically the random effect
variance is the same for all clusters. That is
$$\lambda_{\bullet} = V_{12}^{T}\left( \lambda_{1},...,\lambda_{d} \right)^{T}$$
We return to some specific models below, and show how to fit the ACE and
AE model using this set-up.

One last option in the model-specification is to specify how the
parameters $\lambda_{1},...,\lambda_{d}$ are related. We thus can
specify a matrix $M$ of dimension $p \times d$ such that
$$\left( \lambda_{1},...,\lambda_{d} \right)^{T} = M\theta$$ where
$\theta$ is d-dimensional. If $M$ is diagonal we have no restrictions on
parameters.

This parametrization is obtained with the var.par=0 option that thus
estimates $\theta$.

The DEFAULT parametrization instead estimates the variances of the
random effecs (var.par=1) via the parameters
$\nu$$$M\nu = \left( \lambda_{1}/\lambda_{\bullet}^{2},...,\lambda_{d}/\lambda_{\bullet}^{2} \right)^{T}$$

The basic modelling assumption is now that given random effects
$Z = \left( Z_{1},...,Z_{d} \right)$ we have independent probabilites
$$logit\left( P\left( Y_{i} = 1|X_{i},Z \right) \right) = exp\left( - V_{12,i}^{T}Z \cdot \Psi^{- 1}\left( \lambda_{\bullet},\lambda_{\bullet},P\left( Y_{i} = 1|X_{i} \right) \right) \right)i = 1,2$$

We fit this model by

- Estimating the marginal parameters via GEE
- Using marginal estimates, estimate dependence parameters

Even though the model not formaly in this formulation allows negative
correlation in practice the paramters can be negative and this reflects
negative correlation. An advanatage is that no numerical integration is
needed.

### The twin-stutter data

We consider the twin-stutter where for pairs of twins that are either
dizygotic or monozygotic we have recorded whether the twins are
stuttering

We here consider MZ and same sex DZ twins.

Looking at the data

``` r
library(mets)
data(twinstut)
twinstut$binstut <- 1*(twinstut$stutter=="yes")
twinsall <- twinstut
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
head(twinstut)
#>    tvparnr zyg stutter    sex age nr binstut
#> 1        1  mz      no female  71  1       0
#> 2        1  mz      no female  71  2       0
#> 3        2  dz      no female  71  1       0
#> 8        5  mz      no female  71  1       0
#> 9        5  mz      no female  71  2       0
#> 11       7  dz      no   male  71  1       0
twinstut <- subset(twinstut,tvparnr < 3000)
```

### Pairwise odds ratio model

We start by fitting an overall dependence OR for both MZ and DZ even
though the dependence is expected to be different across zygosity.

The first step is to fit the marginal model adjusting for marginal
covariates. We here note that there is a rather strong gender effect in
the risk of stuttering.

``` r
margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
summary(margbin)
#> 
#> Call:
#> glm(formula = binstut ~ factor(sex) + age, family = binomial(), 
#>     data = twinstut)
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      0.03033    2.16166   0.014   0.9888    
#> factor(sex)male  0.80288    0.20272   3.961 7.48e-05 ***
#> age             -0.05471    0.03295  -1.660   0.0968 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 940.43  on 2650  degrees of freedom
#> Residual deviance: 921.01  on 2648  degrees of freedom
#> AIC: 927.01
#> 
#> Number of Fisher Scoring iterations: 6
```

Now estimating the OR parameter. We see a strong dependence with an OR
at around 8 that is clearly significant.

``` r
bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
                       clusters=twinstut$tvparnr,detail=0)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>                theta       se
#> dependence1 1.859712 0.481585
#> 
#> $or
#>             Estimate Std.Err   2.5% 97.5% P-value
#> dependence1    6.422   3.093 0.3603 12.48 0.03785
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Now, and more interestingly, we consider an OR that depends on zygosity
and note that MZ have a much larger OR than DZ twins. This type of trait
is somewhat complicated to interpret, but clearly, one option is that
that there is a genetic effect, alternatively there might be a stronger
environmental effect for MZ twins.

``` r
# design for OR dependence 
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
bin <- binomial.twostage(margbin,data=twinstut,var.link=1,
                          clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bin)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>                   theta        se
#> factor(zyg)dz 0.1624942 1.0050473
#> factor(zyg)mz 3.7460315 0.7259771
#> 
#> $or
#>               Estimate Std.Err    2.5%   97.5% P-value
#> factor(zyg)dz    1.176   1.182  -1.141   3.494  0.3197
#> factor(zyg)mz   42.353  30.747 -17.910 102.616  0.1684
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

We now consider further regression modelling of the OR structure by
considering possible interactions between sex and zygozsity. We see that
MZ has a much higher dependence and that males have a much lower
dependence. We tested for interaction in this model and these were not
significant.

``` r
twinstut$cage <- scale(twinstut$age)
theta.des <- model.matrix( ~-1+factor(zyg)+factor(sex),data=twinstut)
bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
                          clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>                     theta       se
#> factor(zyg)dz    1.238547 1.118919
#> factor(zyg)mz    6.078248 1.536564
#> factor(sex)male -3.233772 1.776972
#> 
#> $or
#>                  Estimate   Std.Err       2.5%     97.5% P-value
#> factor(zyg)dz     3.45059   3.86094   -4.11670   11.0179  0.3715
#> factor(zyg)mz   436.26434 670.34797 -877.59353 1750.1222  0.5152
#> factor(sex)male   0.03941   0.07003   -0.09784    0.1767  0.5736
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

### Alternative syntax

We now demonstrate how the models can fitted jointly and with anohter
syntax, that ofcourse just fits the marginal model and subsequently fits
the pairwise OR model.

First noticing as before that MZ twins have a much higher dependence.

``` r
 # refers to zygosity of first subject in eash pair : zyg1
 # could also use zyg2 (since zyg2=zyg1 within twinpair's)
 out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,
                response="binstut",id="tvparnr",var.link=1,
                theta.formula=~-1+factor(zyg1))
summary(out)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>                    theta        se
#> factor(zyg1)dz 0.1624942 1.0050473
#> factor(zyg1)mz 3.7460315 0.7259771
#> 
#> $or
#>                Estimate Std.Err    2.5%   97.5% P-value
#> factor(zyg1)dz    1.176   1.182  -1.141   3.494  0.3197
#> factor(zyg1)mz   42.353  30.747 -17.910 102.616  0.1684
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Now considering all data and estimating separate effects for the OR for
opposite sex DZ twins and same sex twins. We here find that os twins are
not markedly different from the same sex DZ twins.

``` r
 # refers to zygosity of first subject in eash pair : zyg1
 # could also use zyg2 (since zyg2=zyg1 within twinpair's))
 
 desfs<-function(x,num1="zyg1",num2="zyg2")
         c(x[num1]=="dz",x[num1]=="mz",x[num1]=="os")*1
     
 margbinall <- glm(binstut~factor(sex)+age,data=twinsall,family=binomial())
 out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
       data=twinsall,response="binstut",id="tvparnr",var.link=1,
       theta.formula=desfs,desnames=c("dz","mz","os"))
 summary(out3)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>        theta        se
#> dz 0.5278527 0.2396796
#> mz 3.4850037 0.1864190
#> os 0.7802940 0.2894394
#> 
#> $or
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> dz    1.695  0.4063  0.8989  2.492 3.016e-05
#> mz   32.623  6.0815 20.7031 44.542 8.128e-08
#> os    2.182  0.6316  0.9442  3.420 5.504e-04
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

### Bivariate Probit model

``` r
library(mets)
data(twinstut)
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
twinstut$binstut <- 1*(twinstut$stutter=="yes")
head(twinstut)
#>    tvparnr zyg stutter    sex age nr binstut
#> 1        1  mz      no female  71  1       0
#> 2        1  mz      no female  71  2       0
#> 3        2  dz      no female  71  1       0
#> 8        5  mz      no female  71  1       0
#> 9        5  mz      no female  71  2       0
#> 11       7  dz      no   male  71  1       0
twinstut <- subset(twinstut,tvparnr < 2000)
```

First testing for same dependence in MZ and DZ that we recommend doing
by comparing the correlations of MZ and DZ twins. Apart from regression
correction in the mean this is an un-structured model, and the useful
concordance and casewise concordance estimates can be reported from this
analysis.

``` r
b1 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="un")
summary(b1)
#> 
#>               Estimate  Std.Err        Z   p-value    
#> (Intercept)   -2.04432  0.10246 -19.9532 < 2.2e-16 ***
#> sexmale        0.41203  0.12392   3.3249 0.0008844 ***
#> atanh(rho) MZ  1.11748  0.43188   2.5875 0.0096683 ** 
#> atanh(rho) DZ  0.23000  0.25889   0.8884 0.3743295    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  Total MZ/DZ Complete pairs MZ/DZ
#>  583/1192    212/354             
#> 
#>                            Estimate 2.5%     97.5%   
#> Tetrachoric correlation MZ  0.80669  0.26456  0.96139
#> Tetrachoric correlation DZ  0.22603 -0.27052  0.62759
#> 
#> MZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.00886  0.00272  0.02840
#> Casewise Concordance  0.43282  0.13214  0.79272
#> Marginal              0.02046  0.01258  0.03312
#> Rel.Recur.Risk       21.15321  2.70681 39.59960
#> log(OR)               4.15336  1.96005  6.34666
#> DZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.00128  0.00013  0.01240
#> Casewise Concordance  0.06251  0.00705  0.38497
#> Marginal              0.02046  0.01258  0.03312
#> Rel.Recur.Risk        3.05481 -3.12067  9.23029
#> log(OR)               1.20535 -1.09090  3.50161
#> 
#>                          Estimate 2.5% 97.5%
#> Broad-sense heritability   1      NaN  NaN
```

### Polygenic modelling

We now turn attention to specific polygenic modelling where special
random effects are used to specify ACE, AE, ADE models and so forth.
This is very easy with the bptwin function. The key parts of the output
are the sizes of the genetic component A and the environmental
component, and we can compare with the results of the unstructed model
above. Also formally we can test if this submodel is acceptable by a
likelihood ratio test.

``` r
b1 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="ace")
summary(b1)
#> 
#>              Estimate   Std.Err       Z p-value   
#> (Intercept)  -4.20958   1.55369 -2.7094 0.00674 **
#> sexmale       0.85990   0.36387  2.3632 0.01812 * 
#> log(var(A))   1.17018   0.99695  1.1738 0.24049   
#> log(var(C)) -23.53038       NaN     NaN     NaN   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  Total MZ/DZ Complete pairs MZ/DZ
#>  583/1192    212/354             
#> 
#>                    Estimate 2.5%     97.5%   
#> A                   0.76318  0.41002  1.11633
#> C                   0.00000  0.00000  0.00000
#> E                   0.23682 -0.11633  0.58998
#> MZ Tetrachoric Cor  0.76318  0.15672  0.95170
#> DZ Tetrachoric Cor  0.38159  0.19280  0.54313
#> 
#> MZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.00768  0.00201  0.02893
#> Casewise Concordance  0.37943  0.09711  0.77658
#> Marginal              0.02025  0.01241  0.03289
#> Rel.Recur.Risk       18.73516 -0.13385 37.60417
#> log(OR)               3.85127  1.58182  6.12072
#> DZ:
#>                      Estimate 2.5%    97.5%  
#> Concordance          0.00230  0.00077 0.00688
#> Casewise Concordance 0.11369  0.05266 0.22841
#> Marginal             0.02025  0.01241 0.03289
#> Rel.Recur.Risk       5.61371  2.19257 9.03484
#> log(OR)              1.92764  1.16374 2.69154
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.76318  0.41002 1.11633
```

``` r
b0 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="ae")
summary(b0)
#> 
#>             Estimate  Std.Err       Z p-value   
#> (Intercept) -4.20958  1.55369 -2.7094 0.00674 **
#> sexmale      0.85990  0.36387  2.3632 0.01812 * 
#> log(var(A))  1.17018  0.99695  1.1738 0.24049   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  Total MZ/DZ Complete pairs MZ/DZ
#>  583/1192    212/354             
#> 
#>                    Estimate 2.5%     97.5%   
#> A                   0.76318  0.41002  1.11633
#> E                   0.23682 -0.11633  0.58998
#> MZ Tetrachoric Cor  0.76318  0.15672  0.95170
#> DZ Tetrachoric Cor  0.38159  0.19280  0.54313
#> 
#> MZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.00768  0.00201  0.02893
#> Casewise Concordance  0.37943  0.09711  0.77658
#> Marginal              0.02025  0.01241  0.03289
#> Rel.Recur.Risk       18.73517 -0.13386 37.60419
#> log(OR)               3.85127  1.58182  6.12072
#> DZ:
#>                      Estimate 2.5%    97.5%  
#> Concordance          0.00230  0.00077 0.00688
#> Casewise Concordance 0.11369  0.05266 0.22841
#> Marginal             0.02025  0.01241 0.03289
#> Rel.Recur.Risk       5.61371  2.19257 9.03484
#> log(OR)              1.92764  1.16374 2.69154
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.76318  0.41002 1.11633
```

### Additive gamma random effects

Fitting first a model with different size random effects for MZ and DZ.
We note that as before in the OR and biprobit model the dependence is
much stronger for MZ twins. We also test if these are the same by
parametrizing the OR model with an intercept. This clearly shows a
significant difference.

``` r
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin <- binomial.twostage(margbin,data=twinstut,model="gamma",
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=1,
     theta.des=theta.des)
summary(bintwin)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link
#> $estimates
#>                     theta        se
#> factor(zyg)dz -2.13348460 1.0821489
#> factor(zyg)mz -0.07740026 0.5790839
#> 
#> $vargam
#>               Estimate Std.Err    2.5%  97.5% P-value
#> factor(zyg)dz   0.1184  0.1282 -0.1327 0.3696 0.35544
#> factor(zyg)mz   0.9255  0.5360 -0.1249 1.9760 0.08419
#> 
#> $type
#> [1] "gamma"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

# test for same dependence in MZ and DZ 
theta.des <- model.matrix( ~factor(zyg),data=twinstut)
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin <- binomial.twostage(margbin,data=twinstut,model="gamma",
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=1,
     theta.des=theta.des)
summary(bintwin)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link
#> $estimates
#>                   theta       se
#> (Intercept)   -2.133485 1.082149
#> factor(zyg)mz  2.056084 1.227348
#> 
#> $vargam
#>               Estimate Std.Err     2.5%   97.5% P-value
#> (Intercept)     0.1184  0.1282  -0.1327  0.3696  0.3554
#> factor(zyg)mz   7.8153  9.5921 -10.9849 26.6155  0.4152
#> 
#> $type
#> [1] "gamma"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

### Polygenic modelling

First setting up the random effects design for the random effects and
the the relationship between variance parameters. We see that the
genetic random effect has size one for MZ and 0.5 for DZ subjects, that
have shared and non-shared genetic components with variance 0.5 such
that the total genetic variance is the same for all subjects. The shared
environmental effect is the samme for all. Thus two parameters with
these bands.

``` r
out <- twin.polygen.design(twinstut,id="tvparnr",zygname="zyg",zyg="dz",type="ace")
head(cbind(out$des.rv,twinstut$tvparnr),10)
#>    MZ DZ DZns1 DZns2 env   
#> 1   1  0     0     0   1  1
#> 2   1  0     0     0   1  1
#> 3   0  1     1     0   1  2
#> 8   1  0     0     0   1  5
#> 9   1  0     0     0   1  5
#> 11  0  1     1     0   1  7
#> 12  0  1     1     0   1  8
#> 13  0  1     0     1   1  8
#> 15  0  1     1     0   1 10
#> 18  0  1     1     0   1 12
out$pardes
#>      [,1] [,2]
#> [1,]  1.0    0
#> [2,]  0.5    0
#> [3,]  0.5    0
#> [4,]  0.5    0
#> [5,]  0.0    1
```

Now, fitting the ACE model, we see that the variance of the genetic,
component, is 1.5 and the environmental variance is -0.5. Thus
suggesting that the ACE model does not fit the data. When the random
design is given we automatically use the gamma fralty model.

``` r
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin1 <- binomial.twostage(margbin,data=twinstut,
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin1)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                  theta        se
#> dependence1  1.1700713 0.9350041
#> dependence2 -0.2552978 0.5910697
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5% P-value
#> dependence1   1.2791  0.6039  0.09552 2.4626 0.03416
#> dependence2  -0.2791  0.6039 -1.46265 0.9045 0.64397
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1   0.9148  0.5352 -0.1343 1.964 0.08744
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

For this model we estimate the concordance and casewise concordance as
well as the marginal rates of stuttering for females.

``` r
concordanceTwinACE(bintwin1,type="ace")
#> $MZ
#>                      Estimate  Std.Err     2.5%   97.5%   P-value
#> concordance          0.009738 0.004010 0.001878 0.01760 1.517e-02
#> casewise concordance 0.476141 0.190613 0.102547 0.84974 1.249e-02
#> marginal             0.020452 0.005083 0.010489 0.03041 5.733e-05
#> 
#> $DZ
#>                      Estimate  Std.Err       2.5%   97.5%   P-value
#> concordance          0.001301 0.001158 -0.0009681 0.00357 2.611e-01
#> casewise concordance 0.063604 0.057255 -0.0486137 0.17582 2.666e-01
#> marginal             0.020452 0.005083  0.0104894 0.03041 5.733e-05
```

The E component was not consistent with the fit of the data and we now
consider instead the AE model.

``` r
out <- twin.polygen.design(twinstut,id="tvparnr",zygname="zyg",zyg="dz",type="ae")

bintwin <- binomial.twostage(margbin,data=twinstut,
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 theta        se
#> dependence1 0.8485392 0.4941264
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err 2.5% 97.5% P-value
#> dependence1        1       0    1     1       0
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1   0.8485  0.4941 -0.1199 1.817 0.08593
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Again, the concordance can be computed:

``` r
concordanceTwinACE(bintwin,type="ae")
#> $MZ
#>                      Estimate  Std.Err     2.5%   97.5%   P-value
#> concordance          0.009236 0.003778 0.001831 0.01664 1.450e-02
#> casewise concordance 0.451605 0.189131 0.080916 0.82229 1.695e-02
#> marginal             0.020452 0.005083 0.010489 0.03041 5.733e-05
#> 
#> $DZ
#>                      Estimate   Std.Err      2.5%    97.5%   P-value
#> concordance          0.001966 0.0007099 0.0005741 0.003357 5.628e-03
#> casewise concordance 0.096106 0.0196561 0.0575803 0.134631 1.012e-06
#> marginal             0.020452 0.0050831 0.0104894 0.030415 5.733e-05
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
#> [1] mets_1.3.9
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.5              knitr_1.51             rlang_1.1.7           
#>  [4] xfun_0.55              textshaping_1.0.4      jsonlite_2.0.0        
#>  [7] listenv_0.10.0         future.apply_1.20.1    lava_1.8.2            
#> [10] htmltools_0.5.9        ragg_1.5.0             sass_0.4.10           
#> [13] rmarkdown_2.30         grid_4.5.2             evaluate_1.0.5        
#> [16] jquerylib_0.1.4        fastmap_1.2.0          numDeriv_2016.8-1.1   
#> [19] yaml_2.3.12            mvtnorm_1.3-3          lifecycle_1.0.5       
#> [22] timereg_2.0.7          compiler_4.5.2         codetools_0.2-20      
#> [25] fs_1.6.6               ucminf_1.2.2           htmlwidgets_1.6.4     
#> [28] Rcpp_1.1.1             future_1.68.0          lattice_0.22-7        
#> [31] systemfonts_1.3.1      digest_0.6.39          R6_2.6.1              
#> [34] parallelly_1.46.1      parallel_4.5.2         splines_4.5.2         
#> [37] Matrix_1.7-4           bslib_0.9.0            tools_4.5.2           
#> [40] RcppArmadillo_15.2.3-1 globals_0.18.0         survival_3.8-3        
#> [43] pkgdown_2.2.0          cachem_1.1.0           desc_1.4.3
```
