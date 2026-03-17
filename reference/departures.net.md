# Departures: netsim Module

This function simulates departure for use in
[netsim](https://epimodel.github.io/EpiModel/reference/netsim.md)
simulations.

## Usage

``` r
departures.net(dat, at)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

## Value

The updated `netsim_dat` main list object.

## See also

[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
