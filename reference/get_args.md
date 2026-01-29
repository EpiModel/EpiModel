# Get Arguments from EpiModel Parameterization Functions

Returns a list of argument names and values for use for parameter
processing functions.

## Usage

``` r
get_args(formal.args, dot.args)
```

## Arguments

- formal.args:

  The output of `formals(sys.function())`.

- dot.args:

  The output of `list(...)`.

## Value

A list of argument names and values.
