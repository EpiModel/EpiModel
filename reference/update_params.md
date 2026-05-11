# Update Model Parameters for Stochastic Network Models

Updates epidemic model parameters originally set with
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
and adds new parameters.

## Usage

``` r
update_params(param, new.param.list)
```

## Arguments

- param:

  Object of class `param.net`, output from function of same name.

- new.param.list:

  Named list of new parameters to add to original parameters.

## Value

An updated list object of class `param.net`, which can be passed to the
EpiModel function
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Details

This function can update any original parameters specified with
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
and add new parameters. This function would be used if the inputs to
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
were a long list of fixed model parameters that needed supplemental
replacements or additions for particular model runs (e.g., changing an
intervention efficacy parameter but leaving all other parameters fixed).

The `new.param.list` object should be a named list object containing
named parameters matching those already in `x` (in which case those
original parameter values will be replaced) or not matching (in which
case new parameters will be added to `param`).

`update_params` modifies a `param.net` object **before** it is passed to
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md);
it is the recommended helper for tweaking parameters outside of a
running simulation. To modify parameters from *inside* a module function
during a simulation, use
[`set_param()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
on the live `dat` object instead. There is no analogous
`update_controls()` or `update_inits()` helper: for a `control.net` or
`init.net` object, edit the list directly (e.g., `ctrl$nsteps <- 1000`)
or rebuild it with a fresh call to the constructor. For *scheduled*
mid-simulation changes to parameters or control settings, see the
`.param.updater.list` argument to
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
and the `.control.updater.list` argument to
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).

## See also

[`set_param()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
for editing parameters inside a module function during a simulation;
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
for the parameter constructor;
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
for running a simulation.

## Examples

``` r
x <- param.net(inf.prob = 0.5, act.rate = 2)
y <- list(inf.prob = 0.75, dx.rate = 0.2)
z <- update_params(x, y)
print(z)
#> Fixed Parameters
#> ---------------------------
#> inf.prob = 0.75
#> act.rate = 2
#> dx.rate = 0.2
```
