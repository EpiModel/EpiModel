# Handle the Logging of Traceback and Dumping of Frames on Error

If `control$.traceback.on.error == TRUE`, this function prints the
traceback of the current simulation to STDIN. This is useful when
`ncores > 1` or in HPC settings. If
`control$.dump.frames.on.error == TRUE`, this function saves a debugging
dump for "postmortem debugging". The dumps are named
"dump\_%Y%m%d\_%H%M%S_s.rda" and stored at the root of the working
directory.

## Usage

``` r
netsim_error_logger(dat, s)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- s:

  The number of the simulation that failed

## Value

Nothing, after logging and dumping frames, the function gives the
control back to the general error handler
