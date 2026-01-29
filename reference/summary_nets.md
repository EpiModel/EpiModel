# Extract Summary Statistics of Networks Used in netsim

This function calls `summary` on each network being simulated in
`netsim`, provided `save.nwstats` and `resimulate.network` are both
`TRUE`. It records the statistics represented by `nwstats.formula` in
`dat$stats$nwstats`, where `dat` is the main `netsim_dat` class object.

## Usage

``` r
summary_nets(dat, at)
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
