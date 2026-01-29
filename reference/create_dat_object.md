# Create a Minimal netsim_dat Main List Object for a Network Model

This helper function populates a `netsim_dat` main list object with the
minimal required elements. All parameters are optional. When none are
given the resulting object is only a shell list of class `netsim_dat`
with the different named elements defined as empty lists.

## Usage

``` r
create_dat_object(
  param = list(),
  init = list(),
  control = list(),
  run = list()
)
```

## Arguments

- param:

  An `EpiModel` object of class
  [`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md).

- init:

  An `EpiModel` object of class
  [`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md).

- control:

  An `EpiModel` object of class
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).

- run:

  A `list` that will contains the objects created by
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
  that are required for between step communication. This list must be
  preserved for restarting models.

## Value

A `netsim_dat` main list object.
