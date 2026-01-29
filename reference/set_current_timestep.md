# Set the Current Timestep

Changes the current timestep in the `netsim_dat` object. Use with
caution. This function exists to work around unforeseen corner cases. In
most situation, `increment_timestep` is preferred.

## Usage

``` r
set_current_timestep(dat, timestep)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- timestep:

  The new value for the timestep.

## Value

The updated `netsim_dat` main list object.

## Mutability

This DOES NOT modify the `netsim_dat` object in place. The result must
be assigned back to `dat` in order to be registered:
`dat <- increment_timestep(dat)`.
