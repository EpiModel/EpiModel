# Increment the Current Timestep

This function adds 1 to the timestep counter stored in the `netsim_dat`
main list object.

## Usage

``` r
increment_timestep(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## Mutability

This DOES NOT modify the `netsim_dat` object in place. The result must
be assigned back to `dat` in order to be registered:
`dat <- increment_timestep(dat)`.
