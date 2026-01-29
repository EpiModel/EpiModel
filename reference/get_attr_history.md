# Extract the Attributes History from Network Simulations

Extract the Attributes History from Network Simulations

## Usage

``` r
get_attr_history(sims)
```

## Arguments

- sims:

  An `EpiModel` object of class `netsim`.

## Value

A list of `data.frame`s, one for each "measure" recorded in the
simulation by the `record_attr_history` function.

## Examples

``` r
if (FALSE) { # \dontrun{

# With `sims` the result of a `netsim` call
get_attr_history(sims)

} # }
```
