# Control Settings for Stochastic Individual Contact Models

Sets the controls for stochastic individual contact models simulated
with [`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
control.icm(type, nsteps, nsims = 1, verbose = FALSE, verbose.int = 0)
```

## Arguments

- type:

  Disease type to be modeled, with the choice of `"SI"` for
  Susceptible-Infected diseases, `"SIR"` for
  Susceptible-Infected-Recovered diseases, and `"SIS"` for
  Susceptible-Infected-Susceptible diseases.

- nsteps:

  Number of time steps to solve the model over. This must be a positive
  integer.

- nsims:

  Number of simulations to run.

- verbose:

  If `TRUE`, print model progress to the console.

- verbose.int:

  Time step interval for printing progress to console, where 0 (the
  default) prints completion status of entire simulation and positive
  integer `x` prints progress after every `x` time steps.

## Value

An `EpiModel` object of class `control.icm`.

## Details

`control.icm` sets the required control settings for any stochastic
individual contact model solved with the
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md)
function. ICM simulations use the built-in SI, SIR, and SIS disease
types only. The `type` argument is required and has no default. For
custom or extension epidemic models, use the network model class via
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
instead.

## See also

Use
[`param.icm()`](https://epimodel.github.io/EpiModel/reference/param.icm.md)
to specify model parameters and
[`init.icm()`](https://epimodel.github.io/EpiModel/reference/init.icm.md)
to specify the initial conditions. Run the parameterized model with
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md).
