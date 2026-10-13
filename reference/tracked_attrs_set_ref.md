# Snapshot Tracked Attributes at Step Start

Captures the current values of all attributes listed in
`control$tracked.attributes` and `control$tracked.attributes.once` along
with their `unique_id`s. This snapshot is stored in
`dat$run$tracked_attrs_ref` and used by `tracked_attrs_record` at the
end of the step to compute deltas. On its first call (when
`dat$run$tracking_attrs` is `FALSE`), records a full initial snapshot
instead for the step `at - 1L`. Attributes that do not exist yet (e.g.,
created by a module on its first call) are left out of the snapshot and
recorded in full at the end of the step in which they appear.

## Usage

``` r
tracked_attrs_set_ref(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## See also

[`tracked_attrs_record()`](https://epimodel.github.io/EpiModel/reference/tracked_attrs_record.md),
[`record_attr_history()`](https://epimodel.github.io/EpiModel/reference/record_attr_history.md)
