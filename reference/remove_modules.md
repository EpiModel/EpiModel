# Remove a Set of Modules From the Module List

Remove a Set of Modules From the Module List

## Usage

``` r
remove_modules(dat, names.to.remove)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- names.to.remove:

  a character vector containing the name of the modules to remove.

## Value

The updated `netsim_dat` main list object.
