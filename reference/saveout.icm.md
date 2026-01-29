# Save icm Data to Output List Format

This function transfers the data from the main `icm_dat` class data
object to the output `out` object at the end of each simulation in
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
saveout.icm(dat, s, out = NULL)
```

## Arguments

- dat:

  Main `icm_dat` class data object passed through `icm` simulations.

- s:

  Current simulation number.

- out:

  Out list passed back in for updating at simulations 2+.

## Value

A list with the following elements:

- **param:** the epidemic parameters passed into the model through
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md),
  with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  [`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md),
  with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model.
