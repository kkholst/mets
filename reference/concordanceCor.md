# Concordance Computes concordance and casewise concordance

Concordance for Twins

## Usage

``` r
concordanceCor(
  object,
  cif1,
  cif2 = NULL,
  messages = TRUE,
  model = NULL,
  coefs = NULL,
  ...
)
```

## Arguments

- object:

  Output from the cor.cif, rr.cif or or.cif function

- cif1:

  Marginal cumulative incidence

- cif2:

  Marginal cumulative incidence of other cause (cause2) if it is
  different from cause1

- messages:

  To print messages

- model:

  Specfifies wich model that is considered if object not given.

- coefs:

  Specfifies dependence parameters if object is not given.

- ...:

  Extra arguments, not used.

## Details

The concordance is the probability that both twins have experienced the
event of interest and is defined as \$\$ cor(t) = P(T_1 \leq t,
\epsilon_1 =1 , T_2 \leq t, \epsilon_2=1) \$\$

Similarly, the casewise concordance is \$\$ casewise(t) =
\frac{cor(t)}{P(T_1 \leq t, \epsilon_1=1) } \$\$ that is the probability
that twin "2" has the event given that twins "1" has.

## References

Estimating twin concordance for bivariate competing risks twin data
Thomas H. Scheike, Klaus K. Holst and Jacob B. Hjelmborg, Statistics in
Medicine 2014, 1193-1204

Estimating Twin Pair Concordance for Age of Onset. Thomas H. Scheike,
Jacob V B Hjelmborg, Klaus K. Holst, 2015 in Behavior genetics
DOI:10.1007/s10519-015-9729-3

## Author

Thomas Scheike
