# Get an Edgelist From the Specified Network

This function outputs an edgelist from the specified network, selecting
the method depending on the stored network type.

## Usage

``` r
get_edgelist(dat, network)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- network:

  Numerical index of the network from which the edgelist should be
  extracted. (May be \> 1 for models with multiple overlapping
  networks.)

## Value

An edgelist in matrix form with two columns. Each column contains the
posit_ids (see `get_posit_ids`) of the nodes in each edge.
