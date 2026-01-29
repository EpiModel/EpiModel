# Copy Vertex Attributes from the `netsim_dat` List to the Network Objects

Copies the vertex attributes stored on the main `attr` list of the
`netsim_dat` object to each of the network objects stored on the
`netsim_dat` object.

## Usage

``` r
copy_datattr_to_nwattr(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## See also

[`get_formula_term_attr()`](http://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md),
[`get_attr_prop()`](http://epimodel.github.io/EpiModel/reference/get_attr_prop.md),
[`auto_update_attr()`](http://epimodel.github.io/EpiModel/reference/auto_update_attr.md),
and
[`copy_nwattr_to_datattr()`](http://epimodel.github.io/EpiModel/reference/copy_nwattr_to_datattr.md).
