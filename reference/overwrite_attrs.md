# Helper to use a `data.frame` to initialize some attributes

Uses `dat$init$init_attr` to overwrite some attributes of the nodes at
initialization

## Usage

``` r
overwrite_attrs(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## Details

If an `init_attr` `data.frame` is present in `dat$init`, use it to
overwrite the attributes it contains. `init_attr` must have a number of
rows equal to the number of nodes in the model as the attributes will be
overwritten one to one, ensuring the correct ordering. `init_attr`
columns MUST have a corresponding attribute already initialized. See
"R/default_attributes.R" for adding new attributes to the model.
`init_attr` is removed from `dat$init` at the end of the function to
free up its memory.
