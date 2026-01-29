# Proportional Table of Vertex Attributes

Calculates the proportional distribution of each vertex attribute
contained in a network.

## Usage

``` r
get_attr_prop(dat, nwterms)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- nwterms:

  Vector of attributes on the network object, usually as output of
  [`get_formula_term_attr()`](http://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md).

## Value

A table containing the proportional distribution of each attribute in
`nwterms`.

## See also

[`get_formula_term_attr()`](http://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md),
[`copy_nwattr_to_datattr()`](http://epimodel.github.io/EpiModel/reference/copy_nwattr_to_datattr.md),
[`auto_update_attr()`](http://epimodel.github.io/EpiModel/reference/auto_update_attr.md).
