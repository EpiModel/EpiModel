# Print Method for ICM Objects

Prints a summary of a stochastic individual contact model object,
including the model type, number of simulations, time steps, model
parameters, and output variable names.

## Usage

``` r
# S3 method for class 'icm'
print(x, ...)
```

## Arguments

- x:

  An object of class `icm`, from
  [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

- ...:

  Additional arguments (currently ignored).

## Details

Given an `icm` object, `print.icm` displays:

- **Simulation summary**: model class, model type (e.g., SI, SIR, SIS),
  number of simulations, number of time steps, and number of groups.

- **Model parameters**: all parameters passed via
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md),
  excluding internal bookkeeping parameters (`groups`, `vital`).

- **Model output**: the names of all epidemic output variables stored in
  the model.
