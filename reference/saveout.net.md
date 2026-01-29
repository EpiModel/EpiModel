# Save netsim Data to Output List Format

This function transfers the data from the main `netsim_dat` object to
the output `out` object at the end of each simulation in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Usage

``` r
saveout.net(dat, s, out = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- s:

  Current simulation number.

- out:

  Out list passed back in for updating at simulations 2+.

## Value

A list with the following elements:

- **param:** the epidemic parameters passed into the model through
  [`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md),
  with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md),
  with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model.

- **stats:** a list containing two sublists, `nwstats` for any network
  statistics saved in the simulation, and `transmat` for the
  transmission matrix saved in the simulation.

- **network:** a list of `networkDynamic` objects, one for each model
  simulation.
