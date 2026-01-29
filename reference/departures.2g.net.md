# Departures: netsim Module

This function simulates departure for use in
[netsim](http://epimodel.github.io/EpiModel/reference/netsim.md)
simulations.

## Usage

``` r
departures.2g.net(dat, at)
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
