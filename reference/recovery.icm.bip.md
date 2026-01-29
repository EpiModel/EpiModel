# Recovery: icm Module

This function simulates recovery from the infected state either to a
distinct recovered state (SIR model type) or back to a susceptible state
(SIS model type), for use in
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
recovery.icm.bip(dat, at)
```

## Arguments

- dat:

  Main `icm_dat` class data object passed through `icm` simulations.

- at:

  Current time step.

## Value

The updated `icm_dat` class main data object.
