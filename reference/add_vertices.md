# Fast Version of network::add.vertices for Edgelist-formatted Network

This function performs a simple operation of updating the edgelist
attribute `n` that tracks the total network size implicit in an edgelist
representation of the network.

## Usage

``` r
add_vertices(el, nv)
```

## Arguments

- el:

  A two-column matrix of current edges (edgelist) with an attribute
  variable `n` containing the total current network size.

- nv:

  A integer equal to the number of nodes to add to the network size at
  the given time step.

## Value

Returns the matrix of current edges, `el`, with the population size
attribute updated based on the number of new vertices specified in `nv`.

## Details

This function is used in `EpiModel` modules to add vertices (nodes) to
the edgelist object to account for entries into the population (e.g.,
births and in-migration).

## Examples

``` r
if (FALSE) { # \dontrun{
library("EpiModel")
nw <- network_initialize(100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

param <- param.net(inf.prob = 0.3)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       tergmLite = TRUE)

# networkLite representation after initialization
dat <- crosscheck.net(x, param, init, control)
dat <- initialize.net(x, param, init, control)

# Check current network size
attributes(dat$el[[1]])$n

# Add 10 vertices
dat$el[[1]] <- add_vertices(dat$el[[1]], 10)

# Check new network size
attributes(dat$el[[1]])$n
} # }
```
