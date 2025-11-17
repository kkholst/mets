# Inverse Probability of Censoring Weights

Internal function. Calculates Inverse Probability of Censoring Weights
(IPCW) and adds them to a data.frame

## Usage

``` r
ipw(
  formula,
  data,
  cluster,
  same.cens = FALSE,
  obs.only = FALSE,
  weight.name = "w",
  trunc.prob = FALSE,
  weight.name2 = "wt",
  indi.weight = "pr",
  cens.model = "aalen",
  pairs = FALSE,
  theta.formula = ~1,
  ...
)
```

## Arguments

- formula:

  Formula specifying the censoring model

- data:

  data frame

- cluster:

  clustering variable

- same.cens:

  For clustered data, should same censoring be assumed (bivariate
  probability calculated as mininum of the marginal probabilities)

- obs.only:

  Return data with uncensored observations only

- weight.name:

  Name of weight variable in the new data.frame

- trunc.prob:

  If TRUE truncation probabilities are also calculated and stored in
  'weight.name2' (based on Clayton-Oakes gamma frailty model)

- weight.name2:

  Name of truncation probabilities

- indi.weight:

  Name of individual censoring weight in the new data.frame

- cens.model:

  Censoring model (default Aalens additive model)

- pairs:

  For paired data (e.g. twins) only the complete pairs are returned
  (With pairs=TRUE)

- theta.formula:

  Model for the dependence parameter in the Clayton-Oakes model
  (truncation only)

- ...:

  Additional arguments to censoring model

## Author

Klaus K. Holst

## Examples

``` r
if (FALSE) { # \dontrun{
data("prt",package="mets")
prtw <- ipw(Surv(time,status==0)~country, data=prt[sample(nrow(prt),5000),],
            cluster="id",weight.name="w")
plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
count <- 0
for (l in unique(prtw$country)) {
    count <- count+1
    prtw <- prtw[order(prtw$time),]
    with(subset(prtw,country==l),
         lines(time,w,col=count,lwd=2))
}
legend("topright",legend=unique(prtw$country),col=1:4,pch=-1,lty=1)
} # }
```
