# Adjustment for the Edges Coefficient with Changing Network Size

Adjusts the edges coefficient in a dynamic network model simulated in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) to
preserve the mean degree of nodes in the network. Requires `at >= 2`.
Maintains the `num(.g2)` epi fields (initialized in
[`sim_nets_t1()`](http://epimodel.github.io/EpiModel/reference/sim_nets_t1.md))
for computing the coefficient adjustment.

## Usage

``` r
edges_correct(dat, at)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

## Value

The updated `netsim_dat` main list object.
