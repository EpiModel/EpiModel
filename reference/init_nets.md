# Network Data and Stats Initialization

This function initializes the network data and stats on the main
`netsim_dat` class data object.

## Usage

``` r
init_nets(dat, x)
```

## Arguments

- dat:

  A main data object of class `netsim_dat` obtained from
  [`create_dat_object()`](http://epimodel.github.io/EpiModel/reference/create_dat_object.md),
  including the `control` argument.

- x:

  Either a fitted network model object of class `netest`, or a list of
  such objects.

## Value

A `netsim_dat` class main data object with network data and stats
initialized.
