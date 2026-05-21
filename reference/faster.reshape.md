# Fast Reshape from Long to Wide Format

Reshapes clustered long-format data to wide format efficiently using
compiled code.

## Usage

``` r
faster.reshape(data, clusters, index.type = FALSE, num = NULL, Rindex = 1)
```

## Arguments

- data:

  a matrix or data.frame to reshape.

- clusters:

  vector of cluster identifiers, or column name in data.

- index.type:

  logical; if TRUE, clusters are already 0-based indices.

- num:

  optional within-cluster numbering variable.

- Rindex:

  if 1, use R (1-based) indexing.

## Value

A wide-format matrix.
