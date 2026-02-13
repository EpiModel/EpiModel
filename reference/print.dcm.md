# Print Method for DCM Objects

Prints a summary of a deterministic compartmental model object,
including the model type, number of runs, time steps, model parameters,
and output variable names.

## Usage

``` r
# S3 method for class 'dcm'
print(x, ...)
```

## Arguments

- x:

  An object of class `dcm`, from
  [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md).

- ...:

  Additional arguments (currently ignored).

## Details

Given a `dcm` object, `print.dcm` displays:

- **Simulation summary**: model class, model type (e.g., SI, SIR, SIS;
  omitted for new/custom models), number of runs, number of time steps,
  and number of groups.

- **Model parameters**: all parameters passed via
  [`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md),
  excluding internal bookkeeping parameters (`groups`, `vital`).

- **Model output**: the names of all epidemic output variables stored in
  the model.
