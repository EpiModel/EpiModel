# Depart Nodes from the netsim_dat Object

Depart Nodes from the netsim_dat Object

## Usage

``` r
depart_nodes(dat, departures)
```

## Arguments

- dat:

  the `netsim_dat` object

- departures:

  the vertex ids of nodes to depart

## Value

the updated `netsim_dat` object with the nodes in `departures` departed

## Details

If `tergmLite` is `FALSE`, the vertex ids `departures` are deactivated
(from the current timestep onward) in each `networkDynamic` stored in
`dat$nw`. If `tergmLite` is `TRUE`, the vertex ids `departures` are
deleted from `dat$el`, `dat$attr`, and `dat$net_attr`.
