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
#>   [1] 51.15051 24.92816 26.46510 52.26570 24.49673 48.73287 50.62500 29.51622
#>   [9] 26.08741 47.34375 22.66165 28.44722 28.01087 35.68808 43.00118 22.61198
#>  [17] 25.72698 35.38619 61.24843 49.03795 60.82419 49.97801 34.33670 41.56035
#>  [25] 39.81226 21.44189 54.89423 22.01149 51.45407 45.06992 41.51404 53.88703
#>  [33] 39.44902 28.43565 28.71728 39.18868 46.72941 30.68254 51.80577 62.40379
#>  [41] 41.01446 42.94289 15.38635 31.93590 60.41651 62.17767 29.89455 46.06695
#>  [49] 38.51219 59.33805 60.78755 43.09923 54.04020 15.20437 46.96594 38.92450
#>  [57] 57.77537 35.05553 23.31573 19.81841 25.72253 44.68851 57.90332 37.16075
#>  [65] 59.68164 21.44140 56.04643 43.15711 49.05869 53.35527 44.33694 21.27178
#>  [73] 50.43239 49.32456 20.37885 28.36782 47.18582 45.81399 40.82687 43.53787
#>  [81] 51.13989 34.54613 62.37561 59.69173 41.29446 54.86499 39.03771 30.41797
#>  [89] 64.04130 43.49662 22.78460 56.32584 41.81394 41.55277 46.53825 31.29447
#>  [97] 58.19170 34.42446 61.96937 18.86591
```
