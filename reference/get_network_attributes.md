# Get Network Attributes from a Network Object

Gets all network attributes except `"mnext"` from its network argument.

## Usage

``` r
get_network_attributes(x)
```

## Arguments

- x:

  An object of class `network` or `networkLite`.

## Value

Returns the named list of network attributes.

## Details

This function is used in `EpiModel` workflows to copy relevant network
attributes from the network object to the `netsim_dat` object when
initializing `netsim` runs.

## Examples

``` r
nw <- network_initialize(100)
get_network_attributes(nw)
#> $bipartite
#> [1] FALSE
#> 
#> $directed
#> [1] FALSE
#> 
#> $hyper
#> [1] FALSE
#> 
#> $loops
#> [1] FALSE
#> 
#> $multiple
#> [1] FALSE
#> 
#> $n
#> [1] 100
#> 
```
