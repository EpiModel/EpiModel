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
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) should be
input into the `init.dcm` function. This function handles initial
conditions for both base model types and original models.

Original models may use the parameter names listed as arguments here, a
new set of names, or a combination of both. With new models, initial
conditions must be input in the same order that the solved derivatives
from the model are output.

## See also

Use
[`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md)
to specify model parameters and
[`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md)
to specify the control settings. Run the parameterized model with
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md).
