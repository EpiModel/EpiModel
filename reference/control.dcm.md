# Control Settings for Deterministic Compartmental Models

Sets the controls for deterministic compartmental models simulated with
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md).

## Usage

``` r
control.dcm(
  type,
  nsteps,
  dt = 1,
  odemethod = "rk4",
  dede = FALSE,
  new.mod = NULL,
  sens.param = TRUE,
  print.mod = FALSE,
  verbose = FALSE,
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

  Number of time steps to solve the model over or vector of times to
  solve the model over. If the number of time steps, then this must be a
  positive integer of length 1.

- dt:

  Time unit for model solutions, with the default of 1. Model solutions
  for fractional time steps may be obtained by setting this to a number
  between 0 and 1.

- odemethod:

  Ordinary differential equation (ODE) integration method, with the
  default of the "Runge-Kutta 4" method (see
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html) for other
  options).

- dede:

  If `TRUE`, use the delayed differential equation solver, which allows
  for time-lagged variables.

- new.mod:

  If not running a base model type, a function with a new model to be
  simulated (see details).

- sens.param:

  If `TRUE`, evaluate arguments in parameters with length greater than 1
  as sensitivity analyses, with one model run per value of the
  parameter. If `FALSE`, one model will be run with parameters of
  arbitrary length (the model may error unless the model function is
  designed to accomodate parameter vectors).

- print.mod:

  If `TRUE`, print the model form to the console.

- verbose:

  If `TRUE`, print model progress to the console.

- ...:

  additional control settings passed to model.

## Value

An `EpiModel` object of class `control.dcm`.

## Details

`control.dcm` sets the required control settings for any deterministic
compartmental models solved with the
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) function.
Controls are required for both base model types and original models. For
all base models, the `type` argument is a necessary parameter and it has
no default.

## New Model Functions

The form of the model function for base models may be displayed with the
`print.mod` argument set to `TRUE`. In this case, the model will not be
run. These model forms may be used as templates to write original model
functions.

These new models may be input and solved with
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) using the
`new.mod` argument, which requires as input a model function.

## See also

Use
[`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md)
to specify model parameters and
[`init.dcm()`](http://epimodel.github.io/EpiModel/reference/init.dcm.md)
to specify the initial conditions. Run the parameterized model with
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md).
