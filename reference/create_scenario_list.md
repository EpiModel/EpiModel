# Make a list of EpiModel scenarios from a data.frame of scenarios

An EpiModel scenario allows one or multiple set of parameters to be
applied to a model a predefined timesteps. They are usually used by a
researcher who wants to model counterfactuals using a pre calibrated
model.

## Usage

``` r
create_scenario_list(scenarios.df)
```

## Arguments

- scenarios.df:

  a `data.frame`

## Value

a list of EpiModel scenarios

## scenarios.df

The `scenarios.df` is a `data.frame` of values to be used as parameters.

It must contain a ".at" column, specifying when the changes should
occur. It requires the "updater" module of EpiModel. *See, vignette*. If
the ".at" value of a row is less than two, the changes will be applied
to the parameter list iteself. The second mandatory column is
".scenario.id". It is used to distinguish the different scenarios. If
multiple rows share the same ".scenario.id", the resulting scenario will
contain one updater per row. This permits modifying parameters at
multiple points in time. (e.g. an intervention limited in time).

The other column names must correspond either to: the name of one
parameter if this parameter is of size 1 or the name of the parameter
with "\_1", "\_2", "*N" with the second part being the position of the
value for a parameter of size \> 1. This means that the parameter names
cannot contain any underscore "*". (e.g "a.rate", "d.rate_1",
"d.rate_2")
