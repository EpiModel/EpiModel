# Initialization: netsim Module

This function initializes the main `netsim_dat` class data object on
which data are stored, simulates the initial state of the networks, and
simulates disease status and other attributes.

## Usage

``` r
initialize.net(x, param, init, control, s)
```

## Arguments

- x:

  If `control$start == 1`, either a fitted network model object of class
  `netest` or a list of such objects. If `control$start > 1`, an object
  of class `netsim`. When multiple networks are used, the node sets
  (including network size and nodal attributes) are assumed to be the
  same for all networks.

- param:

  An `EpiModel` object of class
  [`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md).

- init:

  An `EpiModel` object of class
  [`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md).

- control:

  An `EpiModel` object of class
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).

- s:

  Simulation number, used for restarting dependent simulations.

## Value

A `netsim_dat` class main data object.

## Details

When re-initializing a simulation, the `netsim` object passed to
`initialize.net` must contain the elements `param`, `nwparam`, `epi`,
`coef.form`, and `num.nw`.
