# Adjustment for the Edges Coefficient with Changing Network Size

Adjusts the edges coefficient in a dynamic network model simulated in
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md) to
preserve the mean degree of nodes in the network. Requires `at >= 2`.
Maintains the `num(.g2)` epi fields (initialized in
[`sim_nets_t1()`](https://epimodel.github.io/EpiModel/reference/sim_nets_t1.md))
for computing the coefficient adjustment.

## Usage

``` r
edges_correct(dat, at)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

## Value

The updated `netsim_dat` main list object.

## Details

By default the correction counts every active node. That is correct when
every active node can form a tie, which is the usual case.

It is not correct when a model carries a subpopulation that stays active
but is structurally excluded from the network, for example an age band
past a sexual-cessation age whose target statistics are all zero, so
that `ergm` pins its terms off and no tie incident to those nodes can
form. The population the correction sees then grows while the population
that can hold an edge does not, and the whole of the adjustment lands on
the nodes that can. The `edges.correct.attr` control names a binary
nodal attribute identifying the nodes eligible to form ties, and the
correction then counts only those.

The size of the error is the size of the excluded share. A model whose
excluded band grows to a quarter of the population dilutes mean degree
among the remainder by the same quarter, and none of it appears in the
network diagnostics, because
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) runs
before any node has been excluded.
