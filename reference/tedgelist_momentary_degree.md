# Momentary Degree Distribution of a Timed Edgelist

Momentary Degree Distribution of a Timed Edgelist

## Usage

``` r
tedgelist_momentary_degree(tedgelist, n, window)
```

## Arguments

- tedgelist:

  A timed edgelist, as produced by
  [`networkDynamic::as.data.frame.networkDynamic`](https://rdrr.io/pkg/networkDynamic/man/as.data.frame.networkDynamic.html).

- n:

  Size of the network the edgelist was simulated on.

- window:

  Length-2 numeric of the first and last time step to observe.

## Value

A numeric vector of length `n`, giving the mean number of nodes per time
step with degree 0 through `n - 1` over `window`.
