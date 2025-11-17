# Influence functions or IID decomposition of baseline for recrec/phreg/cifregFG

Influence functions or IID decomposition of baseline for
recrec/phreg/cifregFG

## Usage

``` r
iidBaseline(
  object,
  time = NULL,
  ft = NULL,
  fixbeta = NULL,
  beta.iid = NULL,
  tminus = FALSE,
  ...
)
```

## Arguments

- object:

  phreg/recreg/cifregFG object

- time:

  for baseline IID

- ft:

  function to compute IID of baseline integrated against f(t)

- fixbeta:

  to fix the coefficients

- beta.iid:

  to use these iid of beta

- tminus:

  to get predictions in t-

- ...:

  additional arguments to lower level functions

## Author

Thomas Scheike
