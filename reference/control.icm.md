# Control Settings for Stochastic Individual Contact Models

Sets the controls for stochastic individual contact models simulated
with [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
control.icm(
  type,
  nsteps,
  nsims = 1,
  initialize.FUN = initialize.icm,
  infection.FUN = NULL,
  recovery.FUN = NULL,
  departures.FUN = NULL,
  arrivals.FUN = NULL,
  prevalence.FUN = NULL,
  verbose = FALSE,
  verbose.int = 0,
  skip.check = FALSE,
  ...
)
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

- initialize.FUN:

  Module to initialize the model at the outset, with the default
  function of
  [`initialize.icm()`](http://epimodel.github.io/EpiModel/reference/initialize.icm.md).

- infection.FUN:

  Module to simulate disease infection, with the default function of
  [`infection.icm()`](http://epimodel.github.io/EpiModel/reference/infection.icm.md).

- recovery.FUN:

  Module to simulate disease recovery, with the default function of
  [`recovery.icm()`](http://epimodel.github.io/EpiModel/reference/recovery.icm.md).

- departures.FUN:

  Module to simulate departures or exits, with the default function of
  [`departures.icm()`](http://epimodel.github.io/EpiModel/reference/departures.icm.md).

- arrivals.FUN:

  Module to simulate arrivals or entries, with the default function of
  [`arrivals.icm()`](http://epimodel.github.io/EpiModel/reference/arrivals.icm.md).

- prevalence.FUN:

  Module to calculate disease prevalence at each time step, with the
  default function of
  [`prevalence.icm()`](http://epimodel.github.io/EpiModel/reference/prevalence.icm.md).

- verbose:

  If `TRUE`, print model progress to the console.

- verbose.int:

  Time step interval for printing progress to console, where 0 (the
  default) prints completion status of entire simulation and positive
  integer `x` prints progress after every `x` time steps.

- skip.check:

  If `TRUE`, skips the default error checking for the structure and
  consistency of the parameter values, initial conditions, and control
  settings before running base epidemic models. Setting this to `FALSE`
  is recommended when running models with new modules specified.

- ...:

  Additional control settings passed to model.

## Value

An `EpiModel` object of class `control.icm`.

## Details

`control.icm` sets the required control settings for any stochastic
individual contact model solved with the
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md) function.
Controls are required for both base model types and when passing
original process modules. For all base models, the `type` argument is a
necessary parameter and it has no default.

## New Modules

Base ICM models use a set of module functions that specify how the
individual agents in the population are subjected to infection,
recovery, demographics, and other processes. Core modules are those
listed in the `.FUN` arguments. For each module, there is a default
function used in the simulation. The default infection module, for
example, is contained in the
[`infection.icm()`](http://epimodel.github.io/EpiModel/reference/infection.icm.md)
function.

For original models, one may substitute replacement module functions for
any of the default functions. New modules may be added to the workflow
at each time step by passing a module function via the `...` argument.

## See also

Use
[`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md)
to specify model parameters and
[`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md)
to specify the initial conditions. Run the parameterized model with
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).
