# Progress Print Module for Stochastic Individual Contact Models

This function prints progress from stochastic individual contact models
simulated with
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md) to the
console.

## Usage

``` r
verbose.icm(x, type, s = 1, at = 2)
```

## Arguments

- x:

  If the `type` is "startup", then an object of class `control.icm`;
  otherwise, an object of class `icm_dat`, the main data object in `icm`
  simulations.

- type:

  Progress type, either of "startup" for starting messages before all
  simulations, or "progress" for time step specific messages.

- s:

  Current simulation number, if type is "progress".

- at:

  Current time step, if type is "progress".
