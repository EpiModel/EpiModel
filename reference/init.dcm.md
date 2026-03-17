# Initial Conditions for Deterministic Compartmental Models

Sets the initial conditions for deterministic compartmental models
simulated with `dcm`.

## Usage

``` r
init.dcm(s.num, i.num, r.num, s.num.g2, i.num.g2, r.num.g2, ...)
```

## Arguments

- s.num:

  Number of initial susceptible persons. For two-group models, this is
  the number of initial group 1 susceptible persons.

- i.num:

  Number of initial infected persons. For two-group models, this is the
  number of initial group 1 infected persons.

- r.num:

  Number of initial recovered persons. For two-group models, this is the
  number of initial group 1 recovered persons. This parameter is only
  used for the `SIR` model type.

- s.num.g2:

  Number of initial susceptible persons in group 2. This parameter is
  only used for two-group models.

- i.num.g2:

  Number of initial infected persons in group 2. This parameter is only
  used for two-group models.

- r.num.g2:

  Number of initial recovered persons in group 2. This parameter is only
  used for two-group `SIR` models.

- ...:

  Additional initial conditions passed to model.

## Value

An `EpiModel` object of class `init.dcm`.

## Details

The initial conditions for a model solved with
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md) should
be input into the `init.dcm` function. This function handles initial
conditions for both base model types and original models.

Original models may use the parameter names listed as arguments here, a
new set of names, or a combination of both. With new models, initial
conditions must be input in the same order that the solved derivatives
from the model are output.

## Sensitivity Analyses

Like
[`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md),
initial conditions may be specified as vectors of length greater than
one to run sensitivity analyses over initial conditions. When
`sens.param = TRUE` in
[`control.dcm()`](https://epimodel.github.io/EpiModel/reference/control.dcm.md)
(the default), each element of the vector produces a separate model run.
If both parameters and initial conditions have vector values, all
vectors must have the same length.

## See also

Use
[`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md)
to specify model parameters and
[`control.dcm()`](https://epimodel.github.io/EpiModel/reference/control.dcm.md)
to specify the control settings. Run the parameterized model with
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md).

## Examples

``` r
# SI model initial conditions
init <- init.dcm(s.num = 500, i.num = 1)

# Sensitivity analysis over initial infected count
init <- init.dcm(s.num = 500, i.num = c(1, 5, 25))
```
