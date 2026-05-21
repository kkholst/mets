# Concordance Probability from Twostage Model

Computes concordance probability (joint probability of both subjects
experiencing the event) given dependence parameters and random-effect
variance structures from a twostage model.

Computes concordance probabilities for twin ACE/ADE models from a
binomial twostage object.

Functions for constructing random-effects design matrices for twin and
family models. These designs specify the genetic (A), dominance (D),
common environment (C), and unique environment (E) variance components.

## Usage

``` r
p11_binomial_twostage_RV(
  theta,
  rv1,
  rv2,
  p1,
  p2,
  pardes,
  ags = NULL,
  link = 0,
  i = 1,
  j = 1
)

concordanceTwostage(
  theta,
  p,
  rv1,
  rv2,
  theta.des,
  additive.gamma.sum = NULL,
  link = 0,
  var.par = 0,
  ...
)

concordanceTwinACE(
  object,
  rv1 = NULL,
  rv2 = NULL,
  xmarg = NULL,
  type = "ace",
  ...
)

kendall_ClaytonOakes_twin_ace(parg, parc, K = 10000, test = 0)

kendall.ClaytonOakes.twin.ace(x, y, ...)

kendall_normal_twin_ace(parg, parc, K = 10000)

ascertained_pairs(pairs, data, cr.models, bin = FALSE)

twin.polygen.design(x, ...)

ace_family_design(
  data,
  id = "id",
  member = "type",
  mother = "mother",
  father = "father",
  child = "child",
  child1 = "child",
  type = "ace",
  ...
)

make_pairwise_design(pairs, kinship, type = "ace")
```

## Arguments

- theta:

  dependence parameter vector.

- rv1:

  random-effects design for subject 1.

- rv2:

  random-effects design for subject 2.

- p1:

  marginal probability for subject 1.

- p2:

  marginal probability for subject 2.

- pardes:

  parameter design matrix.

- ags:

  additive gamma sum matrix (optional).

- link:

  link function indicator (0 = identity, 1 = log).

- i:

  index for subject 1.

- j:

  index for subject 2.

- p:

  matrix of marginal probabilities (n x 2).

- theta.des:

  parameter design matrix linking theta to lambda parameters.

- additive.gamma.sum:

  optional matrix for additive gamma sums.

- var.par:

  if 1, parameters are rescaled by sum squared.

- ...:

  additional arguments.

- object:

  a binomial twostage model object.

- xmarg:

  optional covariate values for marginal probabilities.

- type:

  model type: `"ace"`, `"ade"`, `"ae"`, `"de"`, `"dce"`, or `"un"`.

- parg:

  genetic variance parameter (gamma shape for genetic component).

- parc:

  common environment variance parameter (gamma shape for environment).

- K:

  number of simulated twin pairs (multiplied by 2 internally).

- test:

  if 1, prints diagnostic correlations.

- x:

  passed as `parg` (alias wrapper).

- y:

  passed as `parc` (alias wrapper).

- pairs:

  matrix of pair indices (n x 2).

- data:

  a data.frame with twin/family data.

- cr.models:

  formula specifying time and status variables.

- bin:

  logical; if TRUE uses binary (prevalence) ordering rather than time
  ordering.

- id:

  character name of the cluster (pair) identifier column.

- member:

  character name of the family member type column.

- mother:

  value identifying mothers in the member column.

- father:

  value identifying fathers in the member column.

- child:

  value identifying children in the member column.

- child1:

  column name distinguishing first child from second.

- kinship:

  vector of kinship coefficients for each pair.

## Value

A list of concordance tables, one per pair, each containing `pmat` (2x2
probability matrix), `casewise` (casewise concordance), and `marg`
(marginal probabilities).

A list of concordance tables per zygosity group.

A list with components:

- pardes:

  parameter design matrix linking random effects to variance parameters.

- des.rv:

  random-effects design matrix for subjects.

## Details

`twin_polygen_design` creates a polygenic random-effects design for twin
pairs, distinguishing MZ and DZ twins.

`twin.polygen.design` is an alias for `twin_polygen_design`.

`ace_family_design` creates designs for nuclear families (mother,
father, children).

`make_pairwise_design` creates pairwise random-effects designs for
arbitrary kinship structures.

`concordanceTwostage` computes concordance probabilities from a twostage
model.

`concordanceTwinACE` computes concordance from a twin ACE model.

`kendall_ClaytonOakes_twin_ace` and `kendall_normal_twin_ace` compute
Kendall's tau for Clayton-Oakes and normal-frailty twin ACE models
respectively.

`ascertained_pairs` identifies ascertained (affected) pairs in clustered
survival data.

`p11_binomial_twostage_RV` computes the joint probability P(T1\<=t,
T2\<=t) for the additive gamma binary random effects model.

## Author

Klaus K. Holst, Thomas Scheike
