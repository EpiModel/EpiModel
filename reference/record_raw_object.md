# Record an Arbitrary Object During a Simulation

This function records any object during a simulation to allow its
inspection afterward. The records are stored in `dat[["raw.records"]]`
during the simulation, where `dat` is the main `netsim_dat` class
object, and in the `netsim` object under the `raw.records`
[`collections::queue`](https://rdrr.io/pkg/collections/man/queue.html)
object.

## Usage

``` r
record_raw_object(dat, label, object, ..., at = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- label:

  The name to give to the recorded object.

- object:

  The object to be recorded.

- ...:

  Must be empty. It makes `at` a named-only argument, so that calls
  written for the EpiModel 2.6 signature,
  `record_raw_object(dat, at, label, object)`, error instead of silently
  recording the wrong object.

- at:

  The time step at which the recording happens. Defaults to the current
  time step if `NULL`.

## Value

The updated `netsim_dat` main list object.

## Details

See the "Recording Arbitrary Objects (Advanced Debugging)" section of
the "Working with Custom Attributes and Summary Statistics" vignette
([`vignette("attributes-and-summary-statistics", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/attributes-and-summary-statistics.md)).

## Examples

``` r
if (FALSE) { # \dontrun{

dat <- record_raw_object(dat, "a.df", data.frame(x = 2:200))
dat <- record_raw_object(dat, "a.message", "I recorded something")

} # }
```
