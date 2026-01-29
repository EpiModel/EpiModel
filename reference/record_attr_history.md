# Record Attribute History

This function records values specific to a time-step and a group of
nodes. In the records, the `posit_ids` are converted to `unique_ids`
which allows the recording of data for nodes that are no longer in the
network by the end of the run. The records are stored in
`dat[["attr.history"]]` where `dat` is the main `netsim_dat` class
object, and can be accessed from the `netsim` object with
`get_attr_history`.

## Usage

``` r
record_attr_history(dat, at, attribute, posit_ids, values)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  The time where the recording happens.

- attribute:

  The name of the value to record.

- posit_ids:

  A numeric vector of posit_ids to which the measure applies. (see
  `get_posit_ids`).

- values:

  The values to be recorded.

## Value

The updated `netsim_dat` main list object.

## Details

See the "Time-Varying Parameters" section of the "Working With Model
Parameters" vignette.

## Examples

``` r
if (FALSE) { # \dontrun{
# This function must be used inside a custom module
dat <- record_attr_history(dat, at, "attr_1", get_posit_ids(dat), 5)
some_nodes <- get_posit_ids(dat)
some_nodes <- some_nodes[runif(length(some_nodes)) < 0.2]
dat <- record_attr_history(
  dat, at,
  "attr_2",
  some_nodes,
  rnorm(length(some_nodes))
)
} # }
```
