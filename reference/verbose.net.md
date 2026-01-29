# Progress Print Module for Stochastic Network Models

This function prints progress from stochastic network models simulated
with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) to
the console.

## Usage

``` r
verbose.net(x, type, s = 1, at = 2)
```

## Arguments

- x:

  If the `type` is "startup", then an object of class `control.net`;
  otherwise, an object of class `netsim_dat`, the main data object in
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
  simulations.

- type:

  Progress type, either of "startup" for starting messages before all
  simulations, or "progress" for time step specific messages.

- s:

  Current simulation number, if type is "progress".

- at:

  Current time step, if type is "progress".
