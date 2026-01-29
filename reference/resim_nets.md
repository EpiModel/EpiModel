# Resimulate Dynamic Network at Time 2+

This function resimulates the dynamic network in stochastic network
models simulated in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
with dependence between the epidemic and demographic processes and the
network structure.

## Usage

``` r
resim_nets(dat, at)
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
