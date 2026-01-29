# Parameters List for Stochastic Network Models from a Formatted Data Frame

Sets the epidemic parameters for stochastic network models with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
using a specially formatted data frame of parameters.

## Usage

``` r
param.net_from_table(long.param.df)
```

## Arguments

- long.param.df:

  A `data.frame` of parameters. See details for the expected format.

## Value

A list object of class `param.net`, which can be passed to
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## long.param.df

It is possible to set input parameters using a specifically formatted
`data.frame` object. The first 3 columns of this `data.frame` must be:

- `param`: The name of the parameter. If this is a non-scalar parameter
  (a vector of length \> 1), end the parameter name with the position on
  the vector (e.g., `"p_1"`, `"p_2"`, ...).

- `value`: the value for the parameter (or the value of the parameter in
  the Nth position if non-scalar).

- `type`: a character string containing either `"numeric"`, `"logical"`,
  or `"character"` to define the parameter object class.

In addition to these 3 columns, the `data.frame` can contain any number
of other columns, such as `details` or `source` columns to document
parameter meta-data. However, these extra columns will not be used by
EpiModel.
