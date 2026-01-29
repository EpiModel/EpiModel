# Recovery: netsim Module

This function simulates recovery from the infected state either to a
distinct recovered state (SIR model type) or back to a susceptible state
(SIS model type), for use in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Usage

``` r
recovery.net(dat, at)
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
