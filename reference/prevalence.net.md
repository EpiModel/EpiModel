# Get Epidemic Output from netsim Model

Provides all active model state sizes from the network at the specified
time step, output to a list of vectors.

## Usage

``` r
prevalence.net(dat, at)
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

## Details

This network utility is used during the
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
simulation process to efficiently query the current size of each state
or compartment in the model at any given timestep. For a two-group
network, the current state size for each group and overall is provided.
