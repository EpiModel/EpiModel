# Progress Print Module for Deterministic Compartmental Models

This function prints progress from deterministic compartmental models
simulated with
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) to the
console.

## Usage

``` r
verbose.dcm(x, type, s = 1)
```

## Arguments

- x:

  If the `type` is "startup", then an object of class `control.dcm`,
  otherwise the main `df` object in `dcm` runs.

- type:

  Progress type, either of "startup" for starting messages before all
  runs, or "progress" for time step specific messages.

- s:

  Current run number, if type is "progress".
