# Computes mediation weights

Computes mediation weights for either binary or multinomial mediators.
The important part is that the influence functions can be obtained to
compute standard errors.

## Usage

``` r
medweight(
  fit,
  data = data,
  var = NULL,
  name.weight = "weights",
  id.name = "id",
  ...
)
```

## Arguments

- fit:

  either glm-binomial or mlogit (mets package)

- data:

  data frame with data

- var:

  is NULL reads mediator and exposure from formulae in the fit.

- name.weight:

  name of weights

- id.name:

  name of id variable, important for SE computations

- ...:

  Additional arguments to

## Author

Thomas Scheike
