# Specify Controls by Network

This utility function allows specification of certain
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
controls to vary by network. The
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
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

## See also

[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
for passing a list of
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
fits, one per layer. The [Multi-Layer
Networks](https://epimodel.github.io/sismid/11_advanced/mod11-Tutorial.html)
chapter of the Network Modeling for Epidemics course materials works
through a two-layer model in full.
