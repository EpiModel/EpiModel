# Specify Controls by Network

This utility function allows specification of certain
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
controls to vary by network. The
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
control arguments currently supporting `multilayer` specifications are
`nwstats.formula`, `set.control.ergm`, `set.control.tergm`, and
`tergmLite.track.duration`.

## Usage

``` r
multilayer(...)
```

## Arguments

- ...:

  control arguments to apply to each network, with the index of the
  network corresponding to the index of the control argument

## Value

an object of class `multilayer` containing the specified control
arguments
