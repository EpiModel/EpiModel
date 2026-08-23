# Nodes Eligible to Form Ties for the Edges Correction

Returns a logical vector marking the nodes the edges correction should
count, or `NULL` when the `edges.correct.attr` control is unset and
every node counts.

## Usage

``` r
edges_correct_eligible(dat, active)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- active:

  The `active` nodal attribute.

## Value

A logical vector the same length as `active`, or `NULL`.
