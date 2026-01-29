# Deterministic Compartmental Model Functions

These functions parameterize the base deterministic compartmental models
solved using the
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) function.

## Usage

``` r
mod_SI_1g_cl(t, t0, parms)

mod_SI_1g_op(t, t0, parms)

mod_SI_2g_cl(t, t0, parms)

mod_SI_2g_op(t, t0, parms)

mod_SIR_1g_cl(t, t0, parms)

mod_SIR_1g_op(t, t0, parms)

mod_SIR_2g_cl(t, t0, parms)

mod_SIR_2g_op(t, t0, parms)

mod_SIS_1g_cl(t, t0, parms)

mod_SIS_1g_op(t, t0, parms)

mod_SIS_2g_cl(t, t0, parms)

mod_SIS_2g_op(t, t0, parms)
```

## Arguments

- t:

  Time vector, passed into model function internally through
  [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) via the
  control settings in
  [`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md).

- t0:

  Initial conditions for model, passed into model function internally
  through [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md)
  via the initial conditions in
  [`init.dcm()`](http://epimodel.github.io/EpiModel/reference/init.dcm.md).

- parms:

  Model parameters, passed into model function internally through
  [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) via the
  parameter settings in
  [`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md).

## Details

This help page shows the names of all the base deterministic
compartmental model functions supported in EpiModel. Base models are
those already programmed interally within the software. The model
functions may be printed to see their internal structure, either
directly on the console or by using the `print.mod` argument in
[`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md).

The naming convention for the models listed here follows the format:
`mod_<disease type>_<number of groups>_<vital dynamics>`. The supported
disease types are SI, SIS, and SIR; the number of groups are 1 or 2; and
the vital dynamic options are closed (fixed population composition) or
open (with arrivals and departures).
