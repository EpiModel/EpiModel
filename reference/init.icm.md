# Initial Conditions for Stochastic Individual Contact Models

Sets the initial conditions for stochastic individual contact models
simulated with `icm`.

## Usage

``` r
init.icm(s.num, i.num, r.num, s.num.g2, i.num.g2, r.num.g2, ...)
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

An `EpiModel` object of class `init.icm`.

## Details

The initial conditions for a model solved with
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md) should be
input into the `init.icm` function. This function handles initial
conditions for both base models and original models using new modules.

## See also

Use
[`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md)
to specify model parameters and
[`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md)
to specify the control settings. Run the parameterized model with
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).
