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
if (FALSE) { # \dontrun{
library("EpiModel")
set.seed(12345)
nw <- network_initialize(100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

param <- param.net(inf.prob = 0.3)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       tergmLite = TRUE)

# Set seed for reproducibility
set.seed(123456)

# networkLite representation structure after initialization
dat <- crosscheck.net(x, param, init, control)
dat <- initialize.net(x, param, init, control)

# Current edges
head(dat$el[[1]], 20)

# Remove nodes 1 and 2
nodes.to.delete <- 1:2
dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes.to.delete)

# Newly permuted edges
head(dat$el[[1]], 20)
} # }
```
