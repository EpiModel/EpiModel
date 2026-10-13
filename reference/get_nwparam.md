# Extract Network Model Parameters

Extracts a list of network model parameters saved in the initialization
module.

## Usage

``` r
get_nwparam(x, network = 1)
```

## Arguments

- x:

  Main `netsim_dat` data object used in `netsim` simulations.

- network:

  Network number, for simulations with multiple networks representing
  the population.

## Details

On a `netsim` object, the `coef.form` element returned holds the
formation coefficients of its first simulation at the last time step,
which the edges correction may have moved away from the fitted ones.
Those of every simulation are in `x$coef.form`.
