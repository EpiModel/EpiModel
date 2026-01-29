# Get Individual Degree from Network or Edgelist

A fast method for querying the current degree of all individuals within
a network.

## Usage

``` r
get_degree(x)
```

## Arguments

- x:

  Either an object of class `network` or `edgelist` generated from a
  network. If `x` is an edgelist, then it must contain an attribute for
  the total network size, `n`.

## Value

A vector of length equal to the total network size, containing the
current degree of each node in the network.

## Details

Individual-level data on the current degree of nodes within a network is
often useful for summary statistics. Given a `network` class object,
`net`, one way to look up the current degree is to get a summary of the
ERGM term, `sociality`, as in: `summary(net ~ sociality(nodes = NULL))`.
But that is computationally inefficient for a number of reasons. This
function provides a fast method for generating the vector of degrees
using a query of the edgelist. It is even faster if the parameter `x` is
already transformed into an edgelist.

## Examples

``` r
nw <- network_initialize(n = 500)

set.seed(1)
fit <- ergm(nw ~ edges, target.stats = 250)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
#> Evaluating log-likelihood at the estimate. 
#> 
sim <- simulate(fit)

# Slow ERGM-based method
ergm.method <- unname(summary(sim ~ sociality(nodes = NULL)))
ergm.method
#>   [1] 0 1 3 0 1 0 4 0 1 0 3 1 2 1 1 3 0 2 1 0 3 2 3 2 1 1 2 0 2 0 0 0 0 2 2 0 2
#>  [38] 2 4 2 1 0 1 1 1 1 0 2 1 0 1 1 1 0 0 2 2 0 0 1 0 1 2 0 3 1 1 2 0 2 0 0 2 0
#>  [75] 1 0 0 3 0 0 1 0 1 0 2 0 1 1 1 0 3 0 0 1 0 0 0 0 2 1 2 1 0 0 1 0 0 0 0 2 0
#> [112] 0 3 2 1 1 1 1 2 1 0 2 2 0 2 0 0 0 0 1 2 0 1 3 0 0 0 0 2 1 2 0 1 1 2 1 0 1
#> [149] 1 0 0 1 0 0 0 0 3 0 5 2 0 0 0 1 2 2 3 2 0 0 0 1 0 1 1 0 1 1 2 1 0 0 1 2 1
#> [186] 1 0 2 1 0 2 2 2 0 1 0 1 0 1 0 1 1 2 1 0 3 2 2 0 2 0 0 1 1 0 3 1 1 0 0 0 2
#> [223] 2 0 1 0 1 1 1 0 2 2 0 0 1 3 1 0 0 0 1 0 1 1 2 2 0 2 1 1 1 0 2 1 2 2 0 1 0
#> [260] 1 0 0 0 1 1 1 2 0 2 2 1 0 1 1 1 2 0 1 1 0 0 0 4 0 0 0 3 0 2 2 0 2 0 1 0 0
#> [297] 4 1 1 2 0 1 2 0 0 2 0 1 3 1 1 1 1 1 3 1 0 1 1 2 1 0 0 3 0 1 0 1 0 1 1 0 2
#> [334] 2 2 2 1 1 1 0 0 1 0 0 1 0 2 0 1 0 0 1 0 0 2 1 1 2 1 0 2 1 0 0 1 0 0 2 1 2
#> [371] 0 1 1 0 2 2 0 2 0 1 1 0 2 1 1 1 1 0 0 0 4 2 1 1 1 1 3 0 1 0 2 3 2 2 2 0 1
#> [408] 2 0 1 0 2 1 1 0 0 2 1 0 1 0 3 3 2 1 0 2 0 1 2 1 0 1 1 0 3 0 1 0 1 1 4 2 1
#> [445] 2 1 0 0 2 1 1 0 0 1 1 1 0 0 2 2 1 2 1 0 2 0 1 0 0 0 4 2 0 0 1 1 2 0 0 1 1
#> [482] 0 0 0 1 2 0 3 0 2 0 1 1 0 2 3 2 1 0 2

# Fast tabulate method with network object
deg.net <- get_degree(sim)
deg.net
#>   [1] 0 1 3 0 1 0 4 0 1 0 3 1 2 1 1 3 0 2 1 0 3 2 3 2 1 1 2 0 2 0 0 0 0 2 2 0 2
#>  [38] 2 4 2 1 0 1 1 1 1 0 2 1 0 1 1 1 0 0 2 2 0 0 1 0 1 2 0 3 1 1 2 0 2 0 0 2 0
#>  [75] 1 0 0 3 0 0 1 0 1 0 2 0 1 1 1 0 3 0 0 1 0 0 0 0 2 1 2 1 0 0 1 0 0 0 0 2 0
#> [112] 0 3 2 1 1 1 1 2 1 0 2 2 0 2 0 0 0 0 1 2 0 1 3 0 0 0 0 2 1 2 0 1 1 2 1 0 1
#> [149] 1 0 0 1 0 0 0 0 3 0 5 2 0 0 0 1 2 2 3 2 0 0 0 1 0 1 1 0 1 1 2 1 0 0 1 2 1
#> [186] 1 0 2 1 0 2 2 2 0 1 0 1 0 1 0 1 1 2 1 0 3 2 2 0 2 0 0 1 1 0 3 1 1 0 0 0 2
#> [223] 2 0 1 0 1 1 1 0 2 2 0 0 1 3 1 0 0 0 1 0 1 1 2 2 0 2 1 1 1 0 2 1 2 2 0 1 0
#> [260] 1 0 0 0 1 1 1 2 0 2 2 1 0 1 1 1 2 0 1 1 0 0 0 4 0 0 0 3 0 2 2 0 2 0 1 0 0
#> [297] 4 1 1 2 0 1 2 0 0 2 0 1 3 1 1 1 1 1 3 1 0 1 1 2 1 0 0 3 0 1 0 1 0 1 1 0 2
#> [334] 2 2 2 1 1 1 0 0 1 0 0 1 0 2 0 1 0 0 1 0 0 2 1 1 2 1 0 2 1 0 0 1 0 0 2 1 2
#> [371] 0 1 1 0 2 2 0 2 0 1 1 0 2 1 1 1 1 0 0 0 4 2 1 1 1 1 3 0 1 0 2 3 2 2 2 0 1
#> [408] 2 0 1 0 2 1 1 0 0 2 1 0 1 0 3 3 2 1 0 2 0 1 2 1 0 1 1 0 3 0 1 0 1 1 4 2 1
#> [445] 2 1 0 0 2 1 1 0 0 1 1 1 0 0 2 2 1 2 1 0 2 0 1 0 0 0 4 2 0 0 1 1 2 0 0 1 1
#> [482] 0 0 0 1 2 0 3 0 2 0 1 1 0 2 3 2 1 0 2

# Even faster if network already transformed into an edgelist
el <- as.edgelist(sim)
deg.el <- get_degree(el)
deg.el
#>   [1] 0 1 3 0 1 0 4 0 1 0 3 1 2 1 1 3 0 2 1 0 3 2 3 2 1 1 2 0 2 0 0 0 0 2 2 0 2
#>  [38] 2 4 2 1 0 1 1 1 1 0 2 1 0 1 1 1 0 0 2 2 0 0 1 0 1 2 0 3 1 1 2 0 2 0 0 2 0
#>  [75] 1 0 0 3 0 0 1 0 1 0 2 0 1 1 1 0 3 0 0 1 0 0 0 0 2 1 2 1 0 0 1 0 0 0 0 2 0
#> [112] 0 3 2 1 1 1 1 2 1 0 2 2 0 2 0 0 0 0 1 2 0 1 3 0 0 0 0 2 1 2 0 1 1 2 1 0 1
#> [149] 1 0 0 1 0 0 0 0 3 0 5 2 0 0 0 1 2 2 3 2 0 0 0 1 0 1 1 0 1 1 2 1 0 0 1 2 1
#> [186] 1 0 2 1 0 2 2 2 0 1 0 1 0 1 0 1 1 2 1 0 3 2 2 0 2 0 0 1 1 0 3 1 1 0 0 0 2
#> [223] 2 0 1 0 1 1 1 0 2 2 0 0 1 3 1 0 0 0 1 0 1 1 2 2 0 2 1 1 1 0 2 1 2 2 0 1 0
#> [260] 1 0 0 0 1 1 1 2 0 2 2 1 0 1 1 1 2 0 1 1 0 0 0 4 0 0 0 3 0 2 2 0 2 0 1 0 0
#> [297] 4 1 1 2 0 1 2 0 0 2 0 1 3 1 1 1 1 1 3 1 0 1 1 2 1 0 0 3 0 1 0 1 0 1 1 0 2
#> [334] 2 2 2 1 1 1 0 0 1 0 0 1 0 2 0 1 0 0 1 0 0 2 1 1 2 1 0 2 1 0 0 1 0 0 2 1 2
#> [371] 0 1 1 0 2 2 0 2 0 1 1 0 2 1 1 1 1 0 0 0 4 2 1 1 1 1 3 0 1 0 2 3 2 2 2 0 1
#> [408] 2 0 1 0 2 1 1 0 0 2 1 0 1 0 3 3 2 1 0 2 0 1 2 1 0 1 1 0 3 0 1 0 1 1 4 2 1
#> [445] 2 1 0 0 2 1 1 0 0 1 1 1 0 0 2 2 1 2 1 0 2 0 1 0 0 0 4 2 0 0 1 1 2 0 0 1 1
#> [482] 0 0 0 1 2 0 3 0 2 0 1 1 0 2 3 2 1 0 2

identical(as.integer(ergm.method), deg.net, deg.el)
#> [1] TRUE
```
