# Initialize Networks Used in netsim

This function initializes the networks used in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).
The initial edge set for a given network is obtained either from
simulating the cross-sectional model (if `edapprox == TRUE`) or from the
`newnetwork` element of the `netest` object (if `edapprox == FALSE`).
Once the initial edge sets are determined, the first time step is
simulated if `resimulate.network == TRUE`, and all time steps are
simulated if `resimulate.network == FALSE`. Initializes the `num(.g2)`
epi fields used in
[`edges_correct()`](http://epimodel.github.io/EpiModel/reference/edges_correct.md)
for computing edge coefficient adjustments.

## Usage

``` r
sim_nets_t1(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.
