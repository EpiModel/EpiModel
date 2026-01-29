# Get Epidemic Output from icm Model

This function provides all active model state sizes from the network at
the specified time step, output to a list of vectors.

## Usage

``` r
prevalence.icm(dat, at)
```

## Arguments

- dat:

  Main `icm_dat` class data object passed through `icm` simulations.

- at:

  Current time step.

## Value

The updated `icm_dat` class main data object.
