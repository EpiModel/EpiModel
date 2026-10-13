# Record Tracked Attribute Changes at Step End

Compares the current attribute values against the reference snapshot
taken by `tracked_attrs_set_ref` and records only the changes via
`record_attr_history`.

Four types of changes are detected:

- **New attributes**: attributes absent from the reference but present
  now (e.g., created by a module on its first call), recorded for all
  nodes.

- **Departures**: nodes present in the reference but absent now
  (recorded as `active = 0L`).

- **Arrivals**: nodes absent from the reference but present now (all
  tracked attribute values recorded).

- **Value changes**: persisting nodes whose attribute value differs from
  the reference, including `NA` transitions.

A tracked attribute that does not exist at the end of a step is an
error.

## Usage

``` r
tracked_attrs_record(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## See also

[`tracked_attrs_set_ref()`](https://epimodel.github.io/EpiModel/reference/tracked_attrs_set_ref.md),
[`record_attr_history()`](https://epimodel.github.io/EpiModel/reference/record_attr_history.md)
