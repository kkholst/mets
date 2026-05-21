# Coarsen Cluster Identifiers

Reduces the number of unique clusters by binning into quantile groups.

## Usage

``` r
coarse_clust(clusters, max.clust = 100)
```

## Arguments

- clusters:

  vector of cluster identifiers.

- max.clust:

  maximum number of coarsened clusters.

## Value

Integer vector of coarsened cluster indices (0-based).
