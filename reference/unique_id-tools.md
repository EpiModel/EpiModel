# Convert Unique Identifiers to/from Positional Identifiers

EpiModel refers to its nodes either by positional identifiers
(`posit_ids`), which describe the position of a node in the `attr`
vector, or by unique identifiers (`unique_ids`), which allow references
to nodes even after they are deactivated.

## Usage

``` r
get_unique_ids(dat, posit_ids = NULL)

get_posit_ids(dat, unique_ids = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- posit_ids:

  A vector of node positional identifiers (default = NULL).

- unique_ids:

  A vector of node unique identifiers (default = NULL).

## Value

A vector of unique or positional identifiers.

## All elements

When `unique_ids` or `posit_ids` is NULL (default) the full list of
positional IDs or unique IDs is returned.

## Deactivated nodes

When providing `unique_ids` of deactivated nodes to `get_posit_ids`,
`NA`s are returned instead and a warning is produced.
