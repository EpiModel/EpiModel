# Cross Checking of Inputs for Stochastic Network Models

This function checks that the estimation object from
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md) and
the three parameter lists from
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md),
[`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md),
and
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
are consistent.

## Usage

``` r
crosscheck.net(x, param, init, control)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md).

- param:

  An `EpiModel` object of class
  [`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md).

- init:

  An `EpiModel` object of class
  [`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md).

- control:

  An `EpiModel` object of class
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).

## Value

This function returns no objects.
