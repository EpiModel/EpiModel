# Function to Trigger the End Horizon

Function to Trigger the End Horizon

## Usage

``` r
trigger_end_horizon(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## Details

This function triggers the end horizon if a control `end.horizon` exists
and its `at` value is equal to the current timestep. The end horizon
consists on the removal of a set of modules from the module list.
