# Initialization: icm Module

This function initializes the main `icm_dat` class data object, and
simulates disease status and other attributes.

## Usage

``` r
initialize.icm(param, init, control)
```

## Arguments

- param:

  An `EpiModel` object of class
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md).

- init:

  An `EpiModel` object of class
  [`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md).

- control:

  An `EpiModel` object of class
  [`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md).

## Value

The updated `icm_dat` class main data object.
