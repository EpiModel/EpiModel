# Get Vertex Attribute on Network Object

Gets a vertex attribute from an object of class `network`. This
functions simplifies the related function in the `network` package.

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

Returns an object of class `network`.

## Details

This function is used in `EpiModel` workflows to query vertex attributes
on an initialized empty network object (see
[`network_initialize()`](http://epimodel.github.io/EpiModel/reference/network_initialize.md)).

## Examples

``` r
nw <- network_initialize(100)
nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
get_vertex_attribute(nw, "age")
#>   [1] 32.40370 45.68135 42.16617 20.89352 44.00478 30.15790 37.50859 31.57665
#>   [9] 53.61829 21.70735 42.76415 27.23673 43.39668 62.62701 42.56360 32.77265
#>  [17] 27.23223 20.73923 48.55989 33.21743 51.42172 26.12451 45.36705 36.85397
#>  [25] 43.16719 34.05502 49.91213 27.21978 36.93182 17.46908 49.97064 64.76065
#>  [33] 42.98299 46.01470 48.26360 26.84640 49.55004 61.89715 49.93418 26.43900
#>  [41] 64.35403 40.35135 46.35231 42.65546 38.95769 17.43891 60.03029 33.55751
#>  [49] 57.70832 42.73381 16.88279 60.96045 61.01202 29.97121 64.41018 60.77834
#>  [57] 39.68208 28.71390 60.11038 27.61939 51.18840 55.33086 56.68488 64.11624
#>  [65] 26.23454 61.32423 60.81915 16.70723 20.65001 29.14788 45.05377 31.81573
#>  [73] 18.00476 23.05251 50.84553 47.36582 55.06137 31.09444 63.12240 17.62437
#>  [81] 25.21570 39.93163 61.20032 59.35233 48.43742 30.92952 48.48934 29.64794
#>  [89] 36.17831 52.39447 27.05252 31.18440 26.55983 23.25098 23.91550 42.54580
#>  [97] 22.54389 62.94820 22.07813 59.19446
```
