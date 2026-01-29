# Simulate a Network for a Specified Number of Time Steps

This function simulates a dynamic network over one or multiple time
steps for TERGMs or one or multiple cross-sectional network panels for
ERGMs, for use in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
modeling. Network statistics are also extracted and saved if
`save.nwstats == TRUE` and `resimulate.network == FALSE`.

## Usage

``` r
simulate_dat(dat, at, network = 1L, nsteps = 1L)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

- network:

  index of the network to simulate

- nsteps:

  number of time steps to simulate

## Value

The updated `netsim_dat` main list object.
