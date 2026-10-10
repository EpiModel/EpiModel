# Specify Controls and Parameters by Network

This utility function allows specification of certain
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
controls and parameters to vary by network. The
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
control arguments currently supporting `multilayer` specifications are
`nwstats.formula`, `set.control.ergm`, `set.control.tergm`, and
`tergmLite.track.duration`. The
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
parameters supporting them are `inf.prob`, `inf.prob.g2`, and
`act.rate`, which the built-in infection modules then apply per layer.

## Usage

``` r
multilayer(...)
```

## Arguments

- ...:

  control arguments or parameter values to apply to each network, with
  the index of the network corresponding to the index of the argument

## Value

an object of class `multilayer` containing the specified control
arguments or parameter values

## See also

[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
for passing a list of
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
fits and
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
layers, one per network. The [Multi-Layer
Networks](https://epimodel.github.io/sismid/11_advanced/mod11-Tutorial.html)
chapter of the Network Modeling for Epidemics course materials works
through a two-layer model in full.

## Examples

``` r
if (FALSE) { # \dontrun{
# A household layer that transmits more per contact than the community layer
param <- param.net(inf.prob = multilayer(0.45, 0.10), act.rate = 1)
} # }
```
