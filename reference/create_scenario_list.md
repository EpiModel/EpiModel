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

## Vector parameters

A column with a position suffix sets that position of the current value
of the parameter and leaves the other positions as they are. For
example, with `prep.start.rate = c(0.0055, 0.0042, 0.0066)` in `param`,
a scenario with only a `prep.start.rate_1` column of `0.5` gives
`prep.start.rate = c(0.5, 0.0042, 0.0066)`, both when it is applied
before the simulation by
[`use_scenario()`](https://epimodel.github.io/EpiModel/reference/use_scenario.md)
and when it is applied by an updater during the simulation. Before
EpiModel 2.7.0, the parameter was replaced by a vector built from the
scenario columns alone (`0.5` here), so every position had to be given.

A column without a suffix (`prep.start.rate`) still replaces the whole
parameter with its single value. Positions past the end of the current
vector extend it when they leave no gap: with a vector of length 3,
`prep.start.rate_4` alone extends it to length 4, while
`prep.start.rate_5` alone is an error. A scenario cannot shorten a
vector. Only plain vectors can be set by position: a parameter given as
a
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
or
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
object, or as a list, is an error. Use a parameter updater that sets the
whole value instead.
