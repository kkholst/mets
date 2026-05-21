# Create Group Contingency Table from Clustered Data

Creates a contingency table by group from paired/clustered data,
optionally combining lower and upper triangles.

## Usage

``` r
grouptable(
  data,
  id,
  group,
  var,
  lower = TRUE,
  labels,
  order,
  group.labels,
  group.order,
  combine = " & ",
  ...
)
```

## Arguments

- data:

  a data.frame.

- id:

  name of the cluster/pair identifier column.

- group:

  name of the grouping variable (e.g., zygosity).

- var:

  name of the outcome variable to tabulate.

- lower:

  logical; if TRUE, fold upper triangle into lower.

- labels:

  optional labels for levels of `var`.

- order:

  optional ordering of factor levels.

- group.labels:

  optional labels for the groups.

- group.order:

  optional ordering of groups.

- combine:

  separator for combining two groups (default `" & "`).

- ...:

  additional arguments.

## Value

A table or list of tables.
