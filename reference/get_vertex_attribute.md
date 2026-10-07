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
#>   [1] 32.09672 35.55512 20.98375 44.18175 22.29018 38.84287 18.06281 31.46861
#>   [9] 25.53123 47.74077 37.42144 56.33404 16.67101 53.02976 37.14762 42.99052
#>  [17] 32.75700 28.52362 23.52543 31.46900 41.05373 26.22438 39.22338 48.47400
#>  [25] 39.44397 58.88281 28.16707 57.63432 21.83437 59.76652 64.12397 57.12142
#>  [33] 50.93017 50.52766 59.88704 61.71993 16.73064 22.88102 31.89303 52.67740
#>  [41] 38.17126 46.11497 27.55140 54.03835 33.29927 22.69909 54.04413 39.94990
#>  [49] 41.35694 23.97973 29.80143 25.87800 61.98865 21.01547 30.68503 26.11177
#>  [57] 63.78822 53.94593 45.04409 27.62149 63.12511 45.66442 47.97936 55.48544
#>  [65] 37.45336 53.40253 36.39014 39.26885 22.95214 56.69345 24.27507 44.68332
#>  [73] 27.63157 29.43049 48.34490 63.67341 61.22204 45.73420 46.44138 55.98481
#>  [81] 22.00051 15.36000 39.25671 24.01475 38.77540 25.53356 59.30463 54.12232
#>  [89] 23.41538 33.77760 17.61609 64.14295 27.65833 16.26121 15.11226 55.43134
#>  [97] 32.71062 24.34368 19.39466 30.70088
```
