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
must be set separately.

Note that this function only supports arriving new nodes; returning to
an active state nodes that were previously active in the network is not
supported.
