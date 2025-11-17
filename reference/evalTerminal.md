# Evaluates piece constant covariates at min(D,t) where D is a terminal event

returns X(min(D,t)) and min(D,t) and their ratio. for censored
observation 0. to use with the IPCW models implemented.

## Usage

``` r
evalTerminal(
  formula,
  data = data,
  death.code = 2,
  time = NULL,
  marks = NULL,
  mark.codes = NULL
)
```

## Arguments

- formula:

  formula with 'Event' outcome and X to evaluate at min(D,t)

- data:

  data frame

- death.code:

  codes for death (terminating event, 2 default)

- time:

  for evaluation

- marks:

  for terminal events to add marks\*I(D \<=t ,epsilon "in" mark.codes)
  to X(min(D,t))

- mark.codes:

  gives death codes for which to add mark value

## Author

Thomas Scheike
