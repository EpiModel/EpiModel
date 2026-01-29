# Arrivals: netsim Module

This function simulates new arrivals into the network for use in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
simulations.

## Usage

``` r
arrivals.net(dat, at)
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

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
