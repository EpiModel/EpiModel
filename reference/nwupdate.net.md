# Dynamic Network Updates

This function handles all calls to the network object contained on the
main `netsim_dat` object handled in `netsim`.

## Usage

``` r
nwupdate.net(dat, at)
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
