# Set Network State During netsim Simulation

This function updates the `netsim_dat` object given a network
representing the current state of the simulation.

## Usage

``` r
set_network(x, ...)

# S3 method for class 'netsim_dat'
set_network(x, network = 1L, nw, ...)
```

## Arguments

- x:

  a `netsim_dat` object

- network:

  the index of the network to set on `x`

- nw:

  the value of the network to set on `x`

## Value

the `netsim_dat` object with the network state updated

## Details

If running `tergmLite` simulation, this function updates
`x$el[[network]]` and (if `tergmLite.track.duration` is `TRUE` for the
network index `network`) the network attributes `"time"` and
`"lasttoggle"` in `x$net_attr[[network]]`. If not running `tergmLite`
simulation, this function updates the `networkDynamic` object stored in
`x$nw[[network]]`. The input `nw` should be of class `networkLite` when
running `tergmLite` simulation, and of class `networkDynamic` when not
running `tergmLite` simulation.
