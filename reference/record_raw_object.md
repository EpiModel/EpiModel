# Record an Arbitrary Object During a Simulation

This function records any object during a simulation to allow its
inspection afterward. The records are stored in `dat[["raw.records"]]`
during the simulation, where `dat` is the main `netsim_dat` class
object, and in the `netsim` object under the `raw.records`
[`collections::queue`](https://rdrr.io/pkg/collections/man/queue.html)
object.

## Usage

``` r
record_raw_object(dat, at, label, object)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  The time where the recording happens.

- label:

  The name to give to the recorded object.

- object:

  The object to be recorded.

## Value

The updated `netsim_dat` main list object.

## Details

See the "Time-Varying Parameters" section of the "Working With Model
Parameters" vignette.

## Examples

``` r
if (FALSE) { # \dontrun{

dat <- record_raw_object(dat, at, "a.df", data.frame(x = 2:200))
dat <- record_raw_object(dat, at, "a.message", "I recorded something")

} # }
```
