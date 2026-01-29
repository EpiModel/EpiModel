# Group Numbers for Two-Group Network

Outputs group numbers given ID numbers for a two-group network.

## Usage

``` r
idgroup(nw, ids = NULL)
```

## Arguments

- nw:

  Object of class `network` or `networkDynamic`.

- ids:

  Vector of ID numbers for which the group number should be returned. If
  `NULL` (default), return all IDs.

## Value

A vector containing the group number for each of the specified nodes.

## Examples

``` r
nw <- network_initialize(n = 10)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 5))
idgroup(nw)
#>  [1] 1 1 1 1 1 2 2 2 2 2
idgroup(nw, ids = c(3, 6))
#> [1] 1 2
```
