# Set Vertex Attribute on Network Object

Sets a vertex attribute on an object of class `network`. This function
simplifies the related function in the `network` package.

## Usage

``` r
set_vertex_attribute(x, attrname, value, v = NULL)
```

## Arguments

- x:

  An object of class network.

- attrname:

  The name of the attribute to set.

- value:

  A vector of values of the attribute to be set.

- v:

  IDs for the vertices whose attributes are to be altered.

## Value

Returns an object of class `network`.

## Details

This function is used in `EpiModel` workflows to set vertex attributes
on an initialized empty network object (see
[`network_initialize()`](http://epimodel.github.io/EpiModel/reference/network_initialize.md).

## Examples

``` r
nw <- network_initialize(100)
nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
nw
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   total edges= 0 
#>     missing edges= 0 
#>     non-missing edges= 0 
#> 
#>  Vertex attribute names: 
#>     age vertex.names 
#> 
#> No edge attributes
```
