# Cumulative Degree of Each Node in a Timed Edgelist

Cumulative Degree of Each Node in a Timed Edgelist

## Usage

``` r
tedgelist_cumulative_degree(tedgelist, n, window, unique.partners = TRUE)
```

## Arguments

- tedgelist:

  A timed edgelist, as produced by
  [`networkDynamic::as.data.frame.networkDynamic`](https://rdrr.io/pkg/networkDynamic/man/as.data.frame.networkDynamic.html).

- n:

  Size of the network the edgelist was simulated on.

- window:

  Length-2 numeric of the first and last time step to observe.

- unique.partners:

  If `TRUE`, count distinct partners rather than distinct partnerships.

## Value

An integer vector of length `n`, giving the number of partners each node
has over `window`.
