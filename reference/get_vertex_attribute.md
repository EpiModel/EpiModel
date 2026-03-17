# Get Vertex Attribute on Network Object

Gets a vertex attribute from an object of class `network`. This function
simplifies the related function in the `network` package.

## Usage

``` r
get_vertex_attribute(x, attrname)
```

## Arguments

- x:

  An object of class network.

- attrname:

  The name of the attribute to get.

## Value

Returns a vector of vertex attribute values for the attribute specified
by `attrname`.

## Details

This function is used in `EpiModel` workflows to query vertex attributes
on an initialized empty network object (see
[`network_initialize()`](https://epimodel.github.io/EpiModel/reference/network_initialize.md)).

## Examples

``` r
nw <- network_initialize(100)
nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
get_vertex_attribute(nw, "age")
#>   [1] 43.18446 52.02268 23.08326 29.82809 23.53834 27.96412 47.06339 34.05321
#>   [9] 15.36450 25.84393 30.15816 43.67931 54.72058 32.81972 41.62817 32.16732
#>  [17] 63.93147 45.06544 43.45082 16.64340 48.60093 61.49712 30.69253 46.86644
#>  [25] 37.92518 41.87167 46.96175 56.69525 34.72267 31.99487 48.67462 57.67431
#>  [33] 38.58060 38.91728 35.65662 21.80262 58.11414 53.04670 19.25438 60.82528
#>  [41] 30.41047 24.19338 40.87094 29.60332 42.06039 64.00766 24.36419 42.57918
#>  [49] 41.74313 64.25065 42.62640 58.70364 64.15349 16.27865 35.78527 19.73172
#>  [57] 38.17145 18.22314 49.50566 61.96927 16.56265 53.55983 51.09118 50.13094
#>  [65] 24.46165 19.93328 22.31336 15.40808 64.01666 45.86197 42.48395 18.86733
#>  [73] 29.25935 33.04471 45.45428 51.72479 49.60836 48.49588 43.60618 41.03531
#>  [81] 26.13409 43.42104 15.96463 56.64862 29.63464 63.94389 34.67233 62.39618
#>  [89] 35.27721 28.80682 30.58655 49.66932 54.86645 61.87829 56.21041 44.51593
#>  [97] 33.14902 43.93295 28.94766 15.70056
```
