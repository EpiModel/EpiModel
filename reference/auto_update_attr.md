# Update Vertex Attributes for Incoming Vertices

Updates the vertex attributes on a network for new nodes incoming into
that network, based on a set of rules for each attribute that the user
specifies in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).

## Usage

``` r
auto_update_attr(dat, newNodes, curr.tab)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- newNodes:

  Vector of nodal IDs for incoming nodes at the current time step.

- curr.tab:

  Current proportional distribution of all vertex attributes.

## Value

The updated `netsim_dat` main list object.

## See also

[`copy_nwattr_to_datattr()`](http://epimodel.github.io/EpiModel/reference/copy_nwattr_to_datattr.md),
[`get_attr_prop()`](http://epimodel.github.io/EpiModel/reference/get_attr_prop.md),
`auto_update_attr()`.
