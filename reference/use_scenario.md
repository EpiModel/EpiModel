# Apply a scenario object to a param.net object

Apply a scenario object to a param.net object

## Usage

``` r
use_scenario(param, scenario)
```

## Arguments

- param:

  Object of class `param.net`, output from function of same name.

- scenario:

  a scenario object usually created from a `data.frame` of scenarios
  using the `create_scenario_list` function. See the vignette
  "network-model-scenarios".

## Value

An updated list object of class `param.net`, which can be passed to the
EpiModel function
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## scenario

A scenario is a list containing an "id" field, the name of the scenario
and a ".param.updater.list" containing a list of updaters that modifies
the parameters of the model at given time steps. If a scenario contains
a parameter not defined in the `param` object, an error will be
produced. See the vignette "model-parameters" for the technical detail
of their implementation.
