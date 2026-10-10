# Fast Version of network::delete.vertices for Edgelist-formatted Network

Given a current two-column matrix of edges and a vector of IDs to delete
from the matrix, this function first removes any rows of the edgelist in
which the IDs are present and then permutes downward the index of IDs on
the edgelist that were numerically larger than the IDs deleted.

## Usage

``` r
delete_vertices(el, vid)
```

## Arguments

- el:

  A two-column matrix of current edges (edgelist) with an attribute
  variable `n` containing the total current network size.

- vid:

  A vector of IDs to delete from the edgelist.

## Value

Returns an updated edgelist object, `el`, with the edges of deleted
vertices removed from the edgelist and the ID numbers of the remaining
edges permuted downward.

## Details

This function is used in `EpiModel` modules to remove vertices (nodes)
from the edgelist object to account for exits from the population (e.g.,
deaths and out-migration).

## Examples

``` r
library("EpiModel")
set.seed(12345)
nw <- network_initialize(100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

param <- param.net(inf.prob = 0.3)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       tergmLite = TRUE, resimulate.network = TRUE)

# Set seed for reproducibility
set.seed(123456)

# Edgelist representation after initialization
crosscheck.net(x, param, init, control)
dat <- initialize.net(x, param, init, control, s = 1)
el <- get_edgelist(dat, network = 1)

# Current edges
head(el, 20)
#>       [,1] [,2]
#>  [1,]    1   24
#>  [2,]    1   57
#>  [3,]    3   67
#>  [4,]    3   95
#>  [5,]    3   97
#>  [6,]    5   10
#>  [7,]    5   15
#>  [8,]    5   59
#>  [9,]    6   78
#> [10,]    8   25
#> [11,]    9   58
#> [12,]    9   88
#> [13,]   10   14
#> [14,]   10   50
#> [15,]   13   25
#> [16,]   14   65
#> [17,]   18   70
#> [18,]   19   69
#> [19,]   19   93
#> [20,]   19   95

# Remove nodes 1 and 2
nodes.to.delete <- 1:2
el <- delete_vertices(el, nodes.to.delete)

# Newly permuted edges
head(el, 20)
#>       [,1] [,2]
#>  [1,]    1   65
#>  [2,]    1   93
#>  [3,]    1   95
#>  [4,]    3    8
#>  [5,]    3   13
#>  [6,]    3   57
#>  [7,]    4   76
#>  [8,]    6   23
#>  [9,]    7   56
#> [10,]    7   86
#> [11,]    8   12
#> [12,]    8   48
#> [13,]   11   23
#> [14,]   12   63
#> [15,]   16   68
#> [16,]   17   67
#> [17,]   17   91
#> [18,]   17   93
#> [19,]   18   66
#> [20,]   20   40
```
