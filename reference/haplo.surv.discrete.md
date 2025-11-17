# Discrete time to event haplo type analysis

Can be used for logistic regression when time variable is "1" for all
id.

## Usage

``` r
haplo.surv.discrete(
  X = NULL,
  y = "y",
  time.name = "time",
  Haplos = NULL,
  id = "id",
  desnames = NULL,
  designfunc = NULL,
  beta = NULL,
  no.opt = FALSE,
  method = "NR",
  stderr = TRUE,
  designMatrix = NULL,
  response = NULL,
  idhap = NULL,
  design.only = FALSE,
  covnames = NULL,
  fam = binomial,
  weights = NULL,
  offsets = NULL,
  idhapweights = NULL,
  ...
)
```

## Arguments

- X:

  design matrix data-frame (sorted after id and time variable) with id
  time response and desnames

- y:

  name of response (binary response with logistic link) from X

- time.name:

  to sort after time for X

- Haplos:

  (data.frame with id, haplo1, haplo2 (haplotypes (h)) and p=P(h\|G))
  haplotypes given as factor.

- id:

  name of id variale from X

- desnames:

  names for design matrix

- designfunc:

  function that computes design given haplotypes h=(h1,h2) x(h)

- beta:

  starting values

- no.opt:

  optimization TRUE/FALSE

- method:

  NR, nlm

- stderr:

  to return only estimate

- designMatrix:

  gives response and designMatrix directly not implemented (mush
  contain: p, id, idhap)

- response:

  gives response and design directly designMatrix not implemented

- idhap:

  name of id-hap variable to specify different haplotypes for different
  id

- design.only:

  to return only design matrices for haplo-type analyses.

- covnames:

  names of covariates to extract from object for regression

- fam:

  family of models, now binomial default and only option

- weights:

  weights following id for GLM

- offsets:

  following id for GLM

- idhapweights:

  weights following id-hap for GLM (WIP)

- ...:

  Additional arguments to lower level funtions lava::NR optimizer or nlm

## Details

Cycle-specific logistic regression of haplo-type effects with known
haplo-type probabilities. Given observed genotype G and unobserved
haplotypes H we here mix out over the possible haplotypes using that
P(H\|G) is provided.

\$\$ S(t\|x,G)) = E( S(t\|x,H) \| G) = \sum\_{h \in G} P(h\|G) S(t\|z,h)
\$\$ so survival can be computed by mixing out over possible h given g.

Survival is based on logistic regression for the discrete hazard
function of the form \$\$ logit(P(T=t\| T \geq t, x,h)) = \alpha_t +
x(h) \beta \$\$ where x(h) is a regression design of x and haplotypes
\\h=(h_1,h_2)\\

Likelihood is maximized and standard errors assumes that P(H\|G) is
known.

The design over the possible haplotypes is constructed by merging X with
Haplos and can be viewed by design.only=TRUE

## Author

Thomas Scheike

## Examples

``` r
## some haplotypes of interest
types <- c("DCGCGCTCACG","DTCCGCTGACG","ITCAGTTGACG","ITCCGCTGAGG")

## some haplotypes frequencies for simulations 
data(haplo)
hapfreqs <- haplo$hapfreqs 

www <-which(hapfreqs$haplotype %in% types)
hapfreqs$freq[www]
#> [1] 0.138387 0.103394 0.048124 0.291273

baseline=hapfreqs$haplotype[9]
baseline
#> [1] "DTGCGCTCGCG"

designftypes <- function(x,sm=0) {# {{{
hap1=x[1]
hap2=x[2]
if (sm==0) y <- 1*( (hap1==types) | (hap2==types))
if (sm==1) y <- 1*(hap1==types) + 1*(hap2==types)
return(y)
}# }}}

tcoef=c(-1.93110204,-0.47531630,-0.04118204,-1.57872602,-0.22176426,-0.13836416,
0.88830288,0.60756224,0.39802821,0.32706859)

ghaplos <- haplo$ghaplos
haploX  <- haplo$haploX

haploX$time <- haploX$times
Xdes <- model.matrix(~factor(time),haploX)
colnames(Xdes) <- paste("X",1:ncol(Xdes),sep="")
X <- dkeep(haploX,~id+y+time)
X <- cbind(X,Xdes)
Haplos <- dkeep(ghaplos,~id+"haplo*"+p)
desnames=paste("X",1:6,sep="")   # six X's related to 6 cycles 
out <- haplo.surv.discrete(X=X,y="y",time.name="time",
         Haplos=Haplos,desnames=desnames,designfunc=designftypes) 
names(out$coef) <- c(desnames,types)
out$coef
#>          X1          X2          X3          X4          X5          X6 
#> -1.82153345 -0.61608261 -0.17143057 -1.27152045 -0.28635976 -0.19349091 
#> DCGCGCTCACG DTCCGCTGACG ITCAGTTGACG ITCCGCTGAGG 
#>  0.79753613  0.65747412  0.06119231  0.31666905 
summary(out)
#>             Estimate Std.Err     2.5%   97.5%   P-value
#> X1          -1.82153  0.1619 -2.13892 -1.5041 2.355e-29
#> X2          -0.61608  0.1895 -0.98748 -0.2447 1.149e-03
#> X3          -0.17143  0.1799 -0.52398  0.1811 3.406e-01
#> X4          -1.27152  0.2631 -1.78719 -0.7559 1.346e-06
#> X5          -0.28636  0.2030 -0.68425  0.1115 1.584e-01
#> X6          -0.19349  0.2134 -0.61184  0.2249 3.647e-01
#> DCGCGCTCACG  0.79754  0.1494  0.50465  1.0904 9.445e-08
#> DTCCGCTGACG  0.65747  0.1621  0.33971  0.9752 5.007e-05
#> ITCAGTTGACG  0.06119  0.2145 -0.35931  0.4817 7.755e-01
#> ITCCGCTGAGG  0.31667  0.1361  0.04989  0.5834 1.999e-02
```
