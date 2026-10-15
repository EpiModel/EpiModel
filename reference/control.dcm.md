# Control Settings for Deterministic Compartmental Models

Sets the controls for deterministic compartmental models simulated with
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md).

## Usage

``` r
control.dcm(
  type,
  nsteps,
  dt = 1,
  odemethod = "lsoda",
  dede = FALSE,
  new.mod = NULL,
  sens.param = TRUE,
  print.mod = FALSE,
  verbose = FALSE,
  ode.args = NULL,
  solver = NULL,
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
  default of `"lsoda"` (see
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html) for other
  options).

- dede:

  If `TRUE`, use the delayed differential equation solver, which allows
  for time-lagged variables.

- new.mod:

  If not running a base model type, a function with a new model to be
  simulated (see details).

- sens.param:

  If `TRUE`, evaluate arguments in parameters or initial conditions with
  length greater than 1 as sensitivity analyses, with one model run per
  value. If `FALSE`, one model will be run with parameters and initial
  conditions of arbitrary length (the model may error unless the model
  function is designed to accommodate vectors).

- print.mod:

  If `TRUE`, print the model form to the console.

- verbose:

  If `TRUE`, print model progress to the console.

- ode.args:

  A named list of further arguments passed to
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html) (or
  [`deSolve::dede`](https://rdrr.io/pkg/deSolve/man/dede.html) with
  `dede = TRUE`), such as the tolerances `rtol` and `atol` or the
  maximum step `hmax`.

- solver:

  An optional function that solves the model in place of
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html), called as
  `solver(y, times, func, parms, control)` (see details).

- ...:

  additional control settings passed to model.

## Value

An `EpiModel` object of class `control.dcm`.

## Details

`control.dcm` sets the required control settings for any deterministic
compartmental models solved with the
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md)
function. Controls are required for both base model types and original
models. For all base models, the `type` argument is a necessary
parameter and it has no default.

## ODE Solver

As of EpiModel 2.5.5, the default ODE solver was changed from `"rk4"`
(fixed-step Runge-Kutta 4) to `"lsoda"` (adaptive step-size). The
`"lsoda"` method automatically adjusts its internal step size, providing
numerical stability across a wide range of parameter values. The
fixed-step `"rk4"` solver could produce numerical instability (negative
compartment sizes, NaN values) for models with high transmission rates
or stiff dynamics. Users who require the previous behavior can set
`odemethod = "rk4"` explicitly.

## New Model Functions

The form of the model function for base models may be displayed with the
`print.mod` argument set to `TRUE`. In this case, the model will not be
run. These model forms may be used as templates to write original model
functions.

These new models may be input and solved with
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md) using
the `new.mod` argument, which requires as input a model function.

## Custom Solvers

By default,
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md) solves
the model with one call to
[`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html) over all time
steps, passing `odemethod` and any `ode.args`. A model whose rates
change abruptly (an intervention that starts on a given day, for
example) is solved more accurately by restarting the solver at each
change. The `solver` argument allows this, or any other solver, for base
and new models alike. It is called with the initial state `y`, the
solution times `times`, the model function `func`, the parameters
`parms`, and the control settings `control`, and must return a matrix or
data frame with time in the first column and one named column per state
variable, as `ode` does. Any further columns it returns, such as derived
outputs, are kept as model outputs. Columns whose names contain `.flow`
are converted from cumulative counts to counts per time step after the
solver returns, as for the default solver.

## See also

Use
[`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md)
to specify model parameters and
[`init.dcm()`](https://epimodel.github.io/EpiModel/reference/init.dcm.md)
to specify the initial conditions. Run the parameterized model with
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md).
