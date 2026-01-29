# Save dcm Data to Output List Format

This function transfers the data from the main `df` object to the output
`out` object at the end of each run in
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md).

## Usage

``` r
saveout.dcm(df, s, param, control, out = NULL)
```

## Arguments

- df:

  Main object in
  [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md)
  simulations.

- s:

  Current run number.

- param:

  Param list set in
  [`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md).

- control:

  Control list set in
  [`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md).

- out:

  Out list passed back in for updating at runs 2+.

## Value

A list with the following elements:

- **param:** the epidemic parameters passed into the model through
  [`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md),
  with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  [`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md),
  with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model.
