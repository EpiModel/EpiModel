# Create the Unique Identifiers for New Nodes

This function is called by
[`append_core_attr()`](http://epimodel.github.io/EpiModel/reference/net-accessor.md)
and appends new `unique_ids` to the created nodes. It also keeps track
of the already used `unique_ids` with the `dat$run$last_unique_id`
variable.

## Usage

``` r
update_unique_ids(dat, n.new)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- n.new:

  The number of new nodes to give `unique_ids` to.

## Value

The updated `netsim_dat` main list object.
