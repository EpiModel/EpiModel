# Proportional Table of Vertex Attributes

Calculates the proportional distribution of each vertex attribute
contained in a network.

## Usage

``` r
get_attr_prop(dat, nwterms, attrs = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- nwterms:

  Vector of attributes on the network object, usually as output of
  [`get_formula_term_attr()`](https://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md).
  If `NULL`, no table is made.

- attrs:

  Optional character vector of attribute names to restrict the tables
  to. If `NULL` (default), every nodal attribute is tabled.

## Value

A list of proportional tables, one for each nodal attribute (or each one
named in `attrs`), other than `active`, `entrTime`, `exitTime`,
`infTime`, `group`, `status`, `na`, and `vertex.names`. Returns `NULL`
if `nwterms` is `NULL`.

## See also

[`get_formula_term_attr()`](https://epimodel.github.io/EpiModel/reference/get_formula_term_attr.md),
[`copy_nwattr_to_datattr()`](https://epimodel.github.io/EpiModel/reference/copy_nwattr_to_datattr.md),
[`auto_update_attr()`](https://epimodel.github.io/EpiModel/reference/auto_update_attr.md).
