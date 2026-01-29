# Update Model Parameters for Stochastic Network Models

Updates epidemic model parameters originally set with
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
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
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Details

This function can update any original parameters specified with
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
and add new parameters. This function would be used if the inputs to
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
were a long list of fixed model parameters that needed supplemental
replacements or additions for particular model runs (e.g., changing an
intervention efficacy parameter but leaving all other parameters fixed).

The `new.param.list` object should be a named list object containing
named parameters matching those already in `x` (in which case those
original parameter values will be replaced) or not matching (in which
case new parameters will be added to `param`).

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
