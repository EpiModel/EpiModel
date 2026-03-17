# Cross Checking of Inputs for Deterministic Compartmental Models

This function checks that the three parameter lists from
[`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md),
[`init.dcm()`](https://epimodel.github.io/EpiModel/reference/init.dcm.md),
and
[`control.dcm()`](https://epimodel.github.io/EpiModel/reference/control.dcm.md)
are consistent.

## Usage

``` r
crosscheck.dcm(param, init, control)
```

## Arguments

- param:

  An `EpiModel` object of class
  [`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md).

- init:

  An `EpiModel` object of class
  [`init.dcm()`](https://epimodel.github.io/EpiModel/reference/init.dcm.md).

- control:

  An `EpiModel` object of class
  [`control.dcm()`](https://epimodel.github.io/EpiModel/reference/control.dcm.md).

## Value

This function returns no objects.
