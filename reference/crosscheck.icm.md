# Cross Checking of Inputs for Stochastic Individual Contact Models

This function checks that the three parameter lists from
[`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md),
[`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md),
and
[`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md)
are consistent.

## Usage

``` r
crosscheck.icm(param, init, control)
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

This function returns no objects.
