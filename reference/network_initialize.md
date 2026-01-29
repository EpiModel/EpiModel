# Initialize Network Object

Initialize an undirected network object for use in EpiModel workflows.

## Usage

``` r
network_initialize(n)
```

## Arguments

- n:

  Network size.

## Value

Returns an object of class `network`.

## Details

This function is used in `EpiModel` workflows to initialize an empty
network object. The network attributes `directed`, `bipartite`, `hyper`,
`loops`, and `multiple` are set to `FALSE`.

## Examples

``` r
nw <- network_initialize(100)
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
#>     vertex.names 
#> 
#> No edge attributes
```
