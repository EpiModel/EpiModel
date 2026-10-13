# Build a networkDynamic Object from Simulation Output

Reconstructs a `networkDynamic` object from a completed `netsim`
simulation, using the cumulative edgelist for edge spells and the
tracked attribute history for vertex activity and time-varying
attributes. The resulting object can be used with `ndtv` for
visualization or with `tsna` for temporal network analysis.

## Usage

``` r
make_networkDynamic(sim, sim_num = 1, network = NULL)
```

## Arguments

- sim:

  An `EpiModel` object of class `netsim`.

- sim_num:

  The simulation number to use (default 1).

- network:

  Optional network number to filter edges. If `NULL` (default), edges
  from all networks are included.

## Value

A `networkDynamic` object with active vertex spells and time-varying
vertex attributes.

## Details

The simulation must be run with `cumulative.edgelist = TRUE`,
`save.cumulative.edgelist = TRUE`, and `truncate.el.cuml = Inf` in
[`control.net`](https://epimodel.github.io/EpiModel/reference/control.net.md),
so that the cumulative edgelist holds every edge, and with at least one
attribute in `tracked.attributes` or `tracked.attributes.once`, so that
the activity of the nodes is recorded. The vertex ids of the
`networkDynamic` object are the `unique_id`s of the nodes. The edges of
all the network layers are merged into one network, with the layer of
each edge in the `network` edge attribute.

## See also

[`get_attr_at()`](https://epimodel.github.io/EpiModel/reference/get_attr_at.md),
[`get_attr_history()`](https://epimodel.github.io/EpiModel/reference/get_attr_history.md)
