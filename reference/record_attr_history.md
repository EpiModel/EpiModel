# Record Attribute History

Record values for a set of nodes at a specific time step. Records are
stored in `dat[["attr.history"]]` as a
[`collections::queue`](https://rdrr.io/pkg/collections/man/queue.html)
and can be accessed from the `netsim` object with `get_attr_history`.

This function is called automatically by the tracked attributes system
when `tracked.attributes` is set in `control.net`, but can also be
called manually inside custom modules for ad-hoc recording.

## Usage

``` r
record_attr_history(
  dat,
  item,
  value,
  ...,
  at = NULL,
  posit_ids = NULL,
  unique_ids = NULL
)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  The name of the attribute to record.

- value:

  The values to be recorded. Must be of length 1 (recycled for all
  nodes) or the same length as the number of nodes identified by
  `posit_ids` or `unique_ids`.

- ...:

  Must be empty. It makes `at`, `posit_ids` and `unique_ids` named-only
  arguments, so that calls written for the EpiModel 2.6 signature,
  `record_attr_history(dat, at, attribute, posit_ids, values)`, error
  instead of silently recording the wrong values.

- at:

  The time step at which the recording happens. Defaults to the current
  time step if `NULL`.

- posit_ids:

  A numeric vector of positional IDs to which the values apply.
  Converted internally to unique IDs. Either `posit_ids` or `unique_ids`
  must be provided.

- unique_ids:

  A numeric vector of unique IDs to which the values apply. When
  provided, skips the positional-to-unique ID conversion.

## Value

The updated `netsim_dat` main list object.

## Details

Exactly one of `posit_ids` or `unique_ids` must be non-`NULL`. When both
are `NULL`, the function errors. When `unique_ids` is provided directly,
the conversion from positional IDs is skipped, which is more efficient
when unique IDs are already known (e.g., from the automatic tracking
system).

See the "Historical Nodal Attributes" section of the "Working with
Custom Attributes and Summary Statistics" vignette
([`vignette("attributes-and-summary-statistics", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/attributes-and-summary-statistics.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
# Manual recording inside a custom module
dat <- record_attr_history(dat, "status", get_attr(dat, "status"),
                           posit_ids = get_posit_ids(dat))

# Record for a subset of nodes using unique IDs directly
some_uids <- c(10, 25, 42)
dat <- record_attr_history(dat, "risk", c(1, 0, 1),
                           unique_ids = some_uids)
} # }
```
