# Extract Event (Jump) Times

Extracts unique event times from survival data, optionally restricting
to concordant pairs or specific causes.

## Usage

``` r
jumptimes(
  time,
  status = TRUE,
  id,
  cause,
  sample,
  sample.all = TRUE,
  strata = NULL,
  num = NULL,
  ...
)
```

## Arguments

- time:

  vector of event/censoring times.

- status:

  vector of status indicators (default TRUE = event).

- id:

  optional cluster identifier for concordant pair times.

- cause:

  optional cause value to select.

- sample:

  optional maximum number of time points to return.

- sample.all:

  logical; if TRUE and sampling, include remaining times.

- strata:

  not used.

- num:

  optional within-cluster numbering variable.

- ...:

  additional arguments.

## Value

Sorted vector of event times.
