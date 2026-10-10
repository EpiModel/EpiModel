# Coerce a list of parameters to a `long.param.df`

Coerce a list of parameters to a `long.param.df`

## Usage

``` r
param.net_to_table(params)
```

## Arguments

- params:

  A list of parameters to be formatted into a `long.param.df`

## Value

A `data.frame` of parameters.

## Details

Parameters given per network layer with
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
or by duration of infection with
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
are not supported and stop with an error. A parameter table records each
element of a vector as its own row but not what the vector means, so
such a value would read back as a plain vector, which the built-in
modules do not accept. Remove these parameters before converting and set
them again after reading the table back.

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
