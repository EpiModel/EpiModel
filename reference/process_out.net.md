# Save a List of netsim Data to Output List Format

This function transfers the data from a list of the main `netsim_dat`
objects to the output `out` object at the end of all simulations in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Usage

``` r
process_out.net(dat_list)
```

## Arguments

- dat_list:

  A list of main `netsim_dat` objects in `netsim` simulations.

## Value

A list of class `netsim` with the following elements:

- **param:** the epidemic parameters passed into the model through
  `param`, with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  `control`, with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model. Outputs for base models always include the size of
  each compartment, as well as flows in, out of, and between
  compartments.

- **stats:** a list containing two sublists, `nwstats` for any network
  statistics saved in the simulation, and `transmat` for the
  transmission matrix saved in the simulation. See
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
  for further details.

- **network:** a list of `networkDynamic` objects, one for each model
  simulation.
