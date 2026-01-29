# Copy Vertex Attributes From Network to `netsim_dat` List

Copies the vertex attributes stored on the network object to the main
`attr` list in the `netsim_dat` data object.

## Usage

``` r
copy_nwattr_to_datattr(dat, nw)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- nw:

  Network from which to copy vertex attributes.

## Value

The updated `netsim_dat` main list object.

## See also

[`get_formula_term_attr()`](http://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md),
[`get_attr_prop()`](http://epimodel.github.io/EpiModel/reference/get_attr_prop.md),
[`auto_update_attr()`](http://epimodel.github.io/EpiModel/reference/auto_update_attr.md),
and
[`copy_datattr_to_nwattr()`](http://epimodel.github.io/EpiModel/reference/copy_datattr_to_nwattr.md).
