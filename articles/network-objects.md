# Working with Network Objects in EpiModel

## Introduction

This vignette covers how to work with network objects, edgelists, and
partnership histories in EpiModel network models with custom extension
modules. It assumes familiarity with setting up and running network
models with
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
and with the extension API. See the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) (NME) course materials
and the [EpiModel Gallery](https://epimodel.github.io/EpiModel-Gallery/)
for background.

For working with nodal attributes and epidemic summary statistics, see
the companion vignette *Working with Custom Attributes and Summary
Statistics in EpiModel*.

## Network Storage Modes

EpiModel supports two storage modes for networks, controlled by the
`tergmLite` parameter in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md):

- **Full mode** (`tergmLite = FALSE`, the default): Networks are stored
  as `networkDynamic` objects, which preserve the complete history of
  edge activations and deactivations. This allows extraction of the full
  dynamic network after simulation. However, `networkDynamic` objects
  consume substantial memory.

- **tergmLite mode** (`tergmLite = TRUE`): Networks are stored as
  lightweight `networkLite` objects containing only the current edgelist
  and nodal attributes. This provides a 20–50x performance improvement
  and much lower memory usage, making it essential for large-scale
  research models. The trade-off is that
  [`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md)
  returns a `networkLite` (a snapshot) rather than a full dynamic
  network history.

Most extension modules work identically under both modes because they
access networks through EpiModel’s accessor functions rather than
manipulating network objects directly.

## Accessing Network Objects

### During Simulation (Inside Modules)

Inside a custom module, use the
[`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md)
and
[`set_network()`](https://epimodel.github.io/EpiModel/reference/set_network.md)
accessors to work with network objects:

``` r
# Get the network for layer 1
nw <- get_network(dat, network = 1)

# After modifying a network, set it back
dat <- set_network(dat, network = 1, nw = nw)
```

These accessors handle the tergmLite vs. full-mode distinction
internally, so your module code works under both storage modes. In full
mode,
[`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md)
returns a `networkDynamic` object; in tergmLite mode, it returns a
`networkLite`.

In practice, you rarely need to access network objects directly.
Instead, use the edgelist accessor functions described below, which also
work correctly under both storage modes.

### After Simulation

After a `netsim` call, extract network objects with
[`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md):

``` r
sim <- netsim(est, param, init, control)

# Extract the network from simulation 1, network layer 1
nw <- get_network(sim, sim = 1, network = 1)

# Collapse to a static cross-section at time step 50 (full mode only)
nw_at_50 <- get_network(sim, sim = 1, collapse = TRUE, at = 50)
```

In full mode,
[`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md)
returns a `networkDynamic` object. In tergmLite mode, it returns a
`networkLite` object representing the final state. The `collapse` and
`at` arguments are only available in full mode.

**Note:** Network objects are only saved in the output when
`save.network` is `TRUE` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
(the default for full mode). In tergmLite mode, there is no network
history to save.

### Transmission Matrix

The transmission matrix records every transmission event during the
simulation:

``` r
transmat <- get_transmat(sim, sim = 1)
```

This returns a `data.frame` with columns including `at` (time step),
`sus` (ID of the newly infected node), `inf` (ID of the infecting node),
`infDur` (duration of infector’s infection), `transProb`, `actRate`, and
`finalProb`. Transmission matrices are saved by default
(`save.transmat = TRUE` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)).

## Current Edgelists

Current edgelists are the set of active partnerships at the present time
step. These are the most commonly used network data structures inside
extension modules.

### Single Network

[`get_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_edgelist.md)
returns the current edgelist for a given network as a two-column matrix
of positional IDs:

``` r
el <- get_edgelist(dat, network = 1)
```

Each row is an active partnership. Column 1 is the positional ID of the
“head” node; column 2 is the “tail” node. This function works
identically under both storage modes.

### Multiple Networks

[`get_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_edgelists_df.md)
combines edgelists from multiple network layers into a single
`data.frame` with a `network` column:

``` r
# All networks
el_all <- get_edgelists_df(dat, networks = NULL)

# Specific networks
el_12 <- get_edgelists_df(dat, networks = c(1, 2))
```

The returned `data.frame` has columns `head`, `tail`, and `network`.

### Discordant Edgelist

The discordant edgelist identifies partnerships where partners have
different values of a status attribute—the key data structure for
modeling transmission. For example, in an SI model, discordant edges are
those where one partner is susceptible and the other is infected:

``` r
disc_el <- get_discordant_edgelist(
  dat,
  status.attr = "status",
  head.status = "i",
  tail.status = "s"
)
```

The returned `data.frame` has columns `head`, `tail`, `head_status`,
`tail_status`, and `network`. Both orderings are captured: if node A
(infected) is partnered with node B (susceptible), the edge appears
regardless of which is the “head” vs “tail” in the underlying network.

See also
[`discord_edgelist()`](https://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
for the original, simpler version of this function used in built-in
models.

## Positional Indexing and Unique IDs

EpiModel uses two ways to reference nodes:

- **By position:** Think of it like a row number in a spreadsheet.
  `dat$attr$active[3]` accesses the third node’s value directly. This is
  the standard way to look up node information and is very fast. In a
  model with 100 nodes, positions range from 1 to 100. When nodes
  depart, they may be dropped from the vectors, freeing their position
  for new arrivals.

- **By `unique_id`:** A globally unique integer attribute assigned to
  each node at creation and never reused. Slower to look up, but allows
  referencing nodes that have already departed. Used by cumulative
  edgelists and attribute histories.

Conversion between the two systems is handled internally by EpiModel.
The
[`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
and
[`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
functions perform the conversion. See
[`help("unique_id-tools", package = "EpiModel")`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
for details.

## Cumulative Edgelist

The cumulative edgelist is a historical record of all edges in a
network, including the time steps when each edge started and stopped.
This allows querying both current and past partnerships—essential for
contact tracing, partnership duration analysis, and reachability
analysis.

### Enabling Cumulative Edgelists

Cumulative edgelist tracking must be explicitly enabled in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md):

``` r
control <- control.net(
  type = "SI",
  nsims = 1,
  nsteps = 100,
  cumulative.edgelist = TRUE,       # Enable tracking during simulation
  save.cumulative.edgelist = TRUE,  # Save in output after simulation
  verbose = FALSE
)
```

Without `cumulative.edgelist = TRUE`, calls to
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
silently do nothing and calls to
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
return an empty `tibble`.

### Updating the Cumulative Edgelist

The cumulative edgelist must be updated at each time step after network
resimulation. In a custom module or at the end of initialization:

``` r
dat <- update_cumulative_edgelist(dat, network = 1, truncate = Inf)
```

The `truncate` argument controls memory usage:

- `truncate = Inf`: Keep the full history of all edges (no removal).
- `truncate = 0` (the default): Keep only currently active edges. Use
  this if you only need to track active edge start times.
- `truncate = N`: Remove edges that ended more than `N` time steps ago.
  This balances historical depth with memory use.

To update all networks in a multi-layer model:

``` r
for (n_network in seq_len(dat$num.nw)) {
  dat <- update_cumulative_edgelist(dat, n_network, truncate = 100)
}
```

### Accessing the Cumulative Edgelist

#### For a Specific Network

``` r
el_cuml <- get_cumulative_edgelist(dat, network = 1)
```

The returned `tibble` has four columns:

1.  `head`: the `unique_id` of the first node.
2.  `tail`: the `unique_id` of the second node.
3.  `start`: the time step when the edge formed.
4.  `stop`: the last time step the edge was active, or `NA` if the edge
    is still active.

An edge with `start = 5` and `stop = 12` existed from steps 5 through
12, inclusive.
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
always returns a `tibble` with these four columns, even if no edges
exist or the cumulative edgelist was not enabled—in those cases, the
`tibble` has zero rows.

#### For Multiple Networks

``` r
el_cumls <- get_cumulative_edgelists_df(dat, networks = NULL)
```

The `networks` argument accepts a vector of network indices or `NULL`
(all networks). The returned `data.frame` adds a `network` column
identifying which network layer each edge belongs to.

### Contact Tracing

[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md)
extracts the partners of specified nodes from the cumulative edgelist:

``` r
partner_list <- get_partners(
    dat,
    index_posit_ids,
    networks = NULL,
    truncate = Inf,
    only.active.nodes = FALSE
)
```

Arguments:

1.  `dat`: the main list object.
2.  `index_posit_ids`: a vector of positional IDs for the nodes of
    interest (the “indexes”).
3.  `networks`: which network layers to search (`NULL` for all).
4.  `truncate`: only include edges that ended within this many steps
    (filters by edge age).
5.  `only.active.nodes`: if `TRUE`, exclude partnerships with inactive
    (departed) nodes.

The output is similar to
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md)
but with columns `index` and `partner` (both containing unique IDs)
instead of `head` and `tail`. Note that indexes are specified by
positional ID but the output uses unique IDs, since partners may include
nodes that have already departed.

### Cumulative Degree

[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md)
counts the number of distinct partners each node has had over the
tracked history:

``` r
cum_degree <- get_cumulative_degree(
    dat,
    index_posit_ids = 1:50,
    networks = NULL,
    truncate = Inf,
    only.active.nodes = FALSE
)
```

This returns a `data.frame` with columns `index_pid` (positional ID) and
`degree` (cumulative partner count). It wraps
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md)
and counts unique partners per index.

## Reachability Analysis

Reachability functions determine which nodes can be connected through
chains of partnerships over a time window. These are useful for outbreak
investigation (forward reachability: who could a node have infected?)
and source tracing (backward reachability: who could have infected a
node?).

These functions operate on cumulative edgelist objects directly—not on
`dat`—and are typically used for post-simulation analysis.

### Forward Reachable Set

``` r
el_cuml <- get_cumulative_edgelist(dat, network = 1)

fwd <- get_forward_reachable(
  el_cuml,
  from_step = 1,
  to_step = 52,
  nodes = c(10, 25, 42),  # NULL for all nodes with edges
  dense_optim = "auto"
)
```

Returns a list with two elements:

- `reached`: a named list where each element contains the set of nodes
  reachable from each index node through chains of partnerships active
  during `[from_step, to_step]`.
- `lengths`: a matrix with one row per node and one column per time step
  (plus an initial column), showing how the reachable set grows over
  time. This allows back-calculating the distance to each reachable
  node.

Nodes are identified by unique ID in the output (named as `node_ID`).
Nodes with no edges during the analysis period are excluded from the
output; their forward reachable set is just themselves (size 1).

### Backward Reachable Set

``` r
bkwd <- get_backward_reachable(
  el_cuml,
  from_step = 1,
  to_step = 52,
  nodes = c(10, 25, 42)
)
```

Same interface and output structure as
[`get_forward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
but follows partnerships backward in time. This answers the question:
which nodes could have reached this node through a chain of partnerships
during the specified period?

### Performance Notes

Both reachability functions use the
[progressr](https://progressr.futureverse.org/) package for progress
reporting. Wrap calls in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html)
to display a progress bar:

``` r
progressr::with_progress({
  fwd <- get_forward_reachable(el_cuml, from_step = 1, to_step = 260)
})
```

These functions are efficient for large networks because they operate on
cumulative edgelists (much smaller than full `networkDynamic` objects).
For sets of more than 5 nodes, they are faster than iterating over
[`tsna::tPath()`](https://rdrr.io/pkg/tsna/man/paths.html). The
`dense_optim` argument controls an adjacency-list optimization that
helps with dense networks; `"auto"` enables it when the number of edges
exceeds the number of nodes.
