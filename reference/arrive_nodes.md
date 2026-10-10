# Arrive New Nodes to the netsim_dat Object

Arrive New Nodes to the netsim_dat Object

## Usage

``` r
arrive_nodes(dat, nArrivals)
```

## Arguments

- dat:

  the `netsim_dat` object

- nArrivals:

  number of new nodes to arrive

## Value

the updated `netsim_dat` object with `nArrivals` new nodes added

## Details

`nArrivals` new nodes are added to the network data stored on the
`netsim_dat` object. If `tergmLite` is `FALSE`, these nodes are
activated from the current timestep onward. Attributes for the new nodes
must be set separately. On each clique layer (see
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)),
the new nodes are then placed under the layer's `arrivals` rule, which
may set the layer's grouping attribute for them and add their edges. A
model with an observed network layer (see
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md))
has a fixed node set, and both this function and
[`depart_nodes()`](https://epimodel.github.io/EpiModel/reference/depart_nodes.md)
stop if asked to change it.

Note that this function only supports arriving new nodes; returning to
an active state nodes that were previously active in the network is not
supported.
