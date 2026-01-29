# Save Transmission Matrix

This function appends the transmission matrix created during
`infection.net` and `infection.2g.net`.

## Usage

``` r
set_transmat(dat, del, at)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- del:

  Discordant edgelist created within
  [`infection.net()`](http://epimodel.github.io/EpiModel/reference/infection.net.md)
  and
  [`infection.2g.net()`](http://epimodel.github.io/EpiModel/reference/infection.2g.net.md).

- at:

  Current time step.

## Value

The updated `netsim_dat` main list object.

## Details

This internal function works within the parent
[`infection.net()`](http://epimodel.github.io/EpiModel/reference/infection.net.md)
functions to save the transmission matrix created at time step `at` to
the main `netsim_dat` class object `dat`.
