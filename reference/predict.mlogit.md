# Predictions from Multinomial Regression

Computes predicted probabilities for the categorical outcome based on a
`mlogit` object. Can return probabilities for all categories or only for
the observed response.

## Usage

``` r
# S3 method for class 'mlogit'
predict(
  object,
  newdata,
  se = TRUE,
  response = TRUE,
  Y = NULL,
  level = 0.95,
  ...
)
```

## Arguments

- object:

  Object of class `"mlogit"`.

- newdata:

  Data frame for prediction. If missing, predictions are made for the
  original data.

- se:

  Logical; if `TRUE`, computes standard errors and confidence intervals.

- response:

  Logical; if `TRUE` (default), returns probabilities only for the
  observed response (if `newdata` contains the response). If `FALSE`,
  returns probabilities for all categories.

- Y:

  Vector of outcome levels to predict probabilities for (optional).

- level:

  Confidence level for intervals (default 0.95).

- ...:

  Additional arguments.

## Value

A matrix or data frame containing:

- pred:

  Predicted probabilities.

- se:

  Standard errors (if `se=TRUE`).

- lower, upper:

  Confidence intervals (if `se=TRUE`).

If `response=FALSE`, columns correspond to the levels of the outcome
factor.

## See also

[`mlogit`](http://kkholst.github.io/mets/reference/mlogit.md)

## Author

Thomas Scheike
