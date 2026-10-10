# Working with Network Objects in EpiModel

## Introduction

This vignette covers how to work with network objects, edgelists,
multi-layer models, clique and observed network layers, and partnership
histories in EpiModel network models with custom extension modules. It
assumes familiarity with setting up and running network models with
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
`network` (the layer on which the transmission occurred), `infDur`
(duration of infector’s infection), `transProb`, `actRate`, and
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
  `get_attr(dat, "active", posit_ids = 3)` accesses the third node’s
  value directly. This is the standard way to look up node information
  and is very fast. In a model with 100 nodes, positions range from 1
  to 100. When nodes depart, they may be dropped from the vectors,
  freeing their position for new arrivals.

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

## Multi-Layer Networks

A multi-layer model has several edge sets over one node set, such as
main and casual sexual partnerships, sexual and needle-sharing
partnerships, or household and community contacts. Each layer estimated
with
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
has its own formation and dissolution model, and infection can be
transmitted across the edges of any layer. The node set, including the
network size and the nodal attributes, is shared by all layers. A layer
may also be a clique layer built with
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
or an observed network wrapped with
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md);
both are described in the sections that follow. The [Multi-Layer
Networks](https://epimodel.github.io/sismid/11_advanced/mod11-Tutorial.html)
chapter of the NME course materials works through a complete two-layer
model; this section summarizes the mechanics.

### Specifying Layers

Each layer is estimated with its own
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
call on the same starting network, so that every layer carries the same
nodes and nodal attributes. The fitted layers are then passed to
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md) as
a list:

``` r

nw <- network_initialize(n = 1000)
nw <- set_vertex_attribute(nw, "race", rep(0:1, each = 500))

est_main <- netest(nw, formation = ~edges + nodematch("race"),
                   target.stats = c(300, 240),
                   coef.diss = dissolution_coefs(~offset(edges), duration = 100))
est_casl <- netest(nw, formation = ~edges + nodematch("race"),
                   target.stats = c(200, 150),
                   coef.diss = dissolution_coefs(~offset(edges), duration = 10))

sim <- netsim(list(est_main, est_casl), param, init, control)
```

The order of the list sets the network index of each layer: here the
main layer is network 1 and the casual layer is network 2. Every
accessor below refers to a layer by this index.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
reads the nodal attributes from the network of the first layer in the
list. Each layer is diagnosed separately with
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md), as
in a single-layer model.

### Per-Layer Controls and Parameters

The
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
function specifies one value per layer, in the order of the layer list.
Four
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
arguments accept it: `nwstats.formula` (the network statistics recorded
for each layer), `set.control.ergm`, `set.control.tergm`, and
`tergmLite.track.duration`. In
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
`inf.prob`, `inf.prob.g2`, and `act.rate` accept it, and the built-in
infection modules then apply each layer’s value to the discordant edges
of that layer. An argument that is not a `multilayer` object applies to
every layer:

``` r

param <- param.net(inf.prob = multilayer(0.2, 0.1),
                   act.rate = multilayer(2, 1))
control <- control.net(type = "SI", nsteps = 100, nsims = 1,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       nwstats.formula = multilayer("formation",
                                                    ~edges + degree(0:3)))
```

Each entry of a `multilayer` parameter may itself vary with the duration
of infection, given as a
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
object: `inf.prob = multilayer(by_infection_duration(c(0.5, 0.2)), 0.1)`
gives the first layer a probability of 0.5 in the first time step of an
infection and 0.2 after.

In a custom infection module,
[`get_param()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
returns the `multilayer` object, a list with one element per layer, and
the `network` column of the discordant edgelist selects the element that
applies to each edge:

``` r

inf.prob <- get_param(dat, "inf.prob")
del <- get_discordant_edgelist(dat, status.attr = "status",
                               head.status = "i", tail.status = "s")
del$transProb <- unlist(inf.prob)[del$network]  # one fixed value per layer
```

### Accessing Layers

The network accessors take the layer index in their `network` argument:

- **Inside a module:** `get_network(dat, network = 2)` and
  `get_edgelist(dat, network = 2)` read one layer, while
  [`get_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_edgelists_df.md)
  and
  [`get_discordant_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_discordant_edgelist.md)
  combine the layers and add a `network` column (see *Current Edgelists*
  above). `dat$num.nw` holds the number of layers.
- **After the simulation:** `get_network(sim, network = 2)` extracts one
  layer, and `get_nwstats(sim, network = 2)`, `print(sim, network = 2)`,
  and `plot(sim, type = "formation", network = 2)` report its network
  statistics.
- **Transmissions:** the transmission matrix records the layer of each
  transmission in its `network` column, so
  `table(get_transmat(sim)$network)` counts the transmissions on each
  layer.

### Dependent Layers

Layers estimated separately are independent: the edges of one layer do
not affect the formation of edges in another. Dependence between layers
enters through a nodal attribute that summarizes one layer and appears
in the formation model of another. Here, the formation of casual
partnerships depends on whether a node has a main partner:

``` r

# nw carries a starting value of deg.main (1 if the node has a main partner)
est_casl <- netest(nw, formation = ~edges + nodematch("race") +
                     nodefactor("deg.main"),
                   target.stats = c(200, 150, 40),
                   coef.diss = dissolution_coefs(~offset(edges), duration = 10))

update_deg_main <- function(dat, at, network) {
  if (network == 1) {
    deg <- get_degree(get_edgelist(dat, network = 1))
    dat <- set_attr(dat, "deg.main", pmin(deg, 1))
  }
  return(dat)
}

control <- control.net(type = "SI", nsteps = 100, nsims = 1,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       dat.updates = update_deg_main)
```

The starting values of the attribute should be consistent with the main
layer, and they must be on the network of the first layer in the list,
from which
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
reads nodal attributes. The attribute must then be kept current as the
layers are resimulated. The `dat.updates` argument of
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
takes a function of `dat`, `at`, and `network` that
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
calls at every time step, including initialization, before the first
layer is resimulated (`network = 0`) and after each layer is resimulated
(`network = 1, 2, ...`). The function above recomputes `deg.main` after
the main layer is redrawn, so that the casual layer is resimulated
against the current main partnerships. When each layer depends on the
other, the function also recomputes the attribute that summarizes the
second layer; the NME chapter covers that case, including how to seed
the starting values with `san()`.

This pattern assumes `tergmLite = TRUE`, under which the resimulation
reads nodal attributes from `dat`. In full mode, each network object
carries its own copy of the nodal attributes, which
[`nwupdate.net()`](https://epimodel.github.io/EpiModel/reference/nwupdate.net.md)
refreshes from `dat` once per time step, so an attribute set in
`dat.updates` does not reach the next layer’s resimulation in the same
time step. Because main partnerships begin and end while casual
partnerships formed under the earlier value of `deg.main` persist, the
cross-layer statistic can drift from its target during the simulation;
recording it through `nwstats.formula` shows by how much.

## Clique Layers

Some contact structures are groups by definition: every pair of people
in a household, classroom, hospital ward, or ship cabin is in contact
for as long as they share it.
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
builds a network layer of this kind from a grouping attribute,
connecting every pair of nodes that share a value so that each group is
a clique.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
accepts the result anywhere in its list of layers, alongside layers
estimated with
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
and the edges of the layer never form or dissolve on their own.

A clique layer is an addition to
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
not a substitute for it. An ERGM remains the model for any network whose
ties depend on nodal and dyadic predictors (degree, mixing by attribute,
clustering), whether those ties turn over quickly, slowly, or not at
all; a long partnership duration in
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
keeps a `netest` layer close to fixed. A clique layer has no
tie-formation process to estimate. An ERGM could reproduce the cliques
only through a `nodematch` term on the grouping attribute targeted at
its maximum, where no finite coefficient exists. See
[`help("netclique", package = "EpiModel")`](https://epimodel.github.io/EpiModel/reference/netclique.md)
for the full rationale.

### Building a Clique Layer

The grouping attribute is a vertex attribute (integer, numeric, or
character, but not a factor) whose values partition the nodes into
groups. Nodes with a missing (`NA`) value belong to no group and are
isolates on the layer. Two helpers build the attribute:

- [`sample_groups()`](https://epimodel.github.io/EpiModel/reference/sample_groups.md)
  builds a population one group at a time from a table of group types,
  such as household compositions by age group, and returns each node’s
  group ID together with the attributes its group type implies.
- [`assign_groups()`](https://epimodel.github.io/EpiModel/reference/assign_groups.md)
  solves the reverse problem. It assigns group IDs to a population whose
  attributes already exist, drawing group sizes from a target
  distribution and keeping dependent members (such as children) in
  groups with an anchor member (such as an adult).

Here, a population of 1,000 is drawn from four household types, and the
same starting network is used for a household clique layer and a
community TERGM layer:

``` r

hh <- sample_groups(1000, c("adult" = 0.25, "adult adult" = 0.35,
                            "adult adult child" = 0.25,
                            "adult adult child child" = 0.15),
                    attr.name = "age")

nw <- network_initialize(n = 1000)
nw <- set_vertex_attribute(nw, "age", hh$age)
nw <- set_vertex_attribute(nw, "hh_id", hh$group)

est_hh <- netclique(nw, group.attr = "hh_id", arrivals = "join")
est_com <- netest(nw, formation = ~edges + nodematch("age"),
                  target.stats = c(400, 300),
                  coef.diss = dissolution_coefs(~offset(edges), duration = 20))

param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
sim <- netsim(list(est_hh, est_com), param, init, control)
```

The household layer is network 1 and the community layer is network 2,
and the
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
parameter gives household contacts the higher per-act transmission
probability. Since
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
reads nodal attributes from the first layer in the list, the grouping
attribute should also be set on the network passed to
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
for the other layers, as it is here; when it is not,
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
copies it from the clique layer.

Printing a `netclique` object shows the number of groups, the edge
count, the mean degree, and the group size distribution, and
`print(est_hh, by = "age")` adds the mean degree by a nodal attribute.
There is no model to diagnose, so
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) does
not accept a `netclique` object.

### Clique Layers During Simulation

A clique layer is stored like any other layer: as an edgelist in
tergmLite mode and as a `networkDynamic` object in full mode. The
accessors in this vignette
([`get_network()`](https://epimodel.github.io/EpiModel/reference/get_network.md),
[`get_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_edgelist.md),
[`get_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_edgelists_df.md),
[`get_discordant_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_discordant_edgelist.md))
read it without special handling, and the built-in infection modules
transmit across its edges as on any other layer. The network
resimulation and the edges correction skip the layer, so its edges
change only when nodes depart, when nodes arrive, and when a module
calls
[`move_to_group()`](https://epimodel.github.io/EpiModel/reference/move_to_group.md).
With the default `nwstats.formula = "formation"`, the network statistic
recorded for the layer is its edge count.

### Arrivals and Departures

Departing nodes are removed from a clique layer together with their
edges, as on every other layer. Arriving nodes are placed on the layer
under the rule chosen with the `arrivals` argument of
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md).
The rule is applied within
[`arrive_nodes()`](https://epimodel.github.io/EpiModel/reference/arrive_nodes.md),
called from the built-in
[`nwupdate.net()`](https://epimodel.github.io/EpiModel/reference/nwupdate.net.md)
module after the arrivals module has created the new nodes and set their
other attributes:

- `"isolate"` (the default): the new node has no edges on the layer and
  its grouping attribute is `NA`. A custom module may place it later
  with
  [`move_to_group()`](https://epimodel.github.io/EpiModel/reference/move_to_group.md).
- `"new"`: each new node starts a group of size one, with a fresh group
  ID.
- `"join"`: each new node joins an existing group and is connected to
  every active member of that group, including other nodes joining the
  same group in the same time step. By default, the group is drawn with
  probability proportional to its current size.

Under `"join"`, an `arrivals.FUN` function replaces the default draw. It
is called as `arrivals.FUN(dat, at, new_ids, network)`, where `new_ids`
holds the positional IDs of the new nodes, and it returns one group ID
per new node (`NA` leaves that node as an isolate). Because it receives
`dat`, it can read any nodal attribute. Here, each newborn joins a
household that already has a child:

``` r

est_hh <- netclique(nw, group.attr = "hh_id", arrivals = "join",
  arrivals.FUN = function(dat, at, new_ids, network) {
    active <- get_attr(dat, "active")
    hh <- get_attr(dat, "hh_id")
    age <- get_attr(dat, "age")
    existing <- active == 1
    existing[new_ids] <- FALSE
    pool <- hh[which(existing & age == "child" & !is.na(hh))]
    pool[sample.int(length(pool), length(new_ids), replace = TRUE)]
  })
```

The pool leaves out `new_ids`. When the arrivals module gives new nodes
a default value of the grouping attribute, such as `0`, a pool taken
from all nodes would contain that value, and every arrival that drew it
would join the other arrivals in a spurious group. Define `arrivals.FUN`
at the top level of a script or in a package: the layer keeps the
function with its enclosing environment, and so does every simulation
run from it.

A custom arrivals module does not need to set the grouping attribute:
the new nodes receive `NA` unless the module gives the attribute a
default value, and the layer’s rule then places them either way. A
module that does set the attribute itself should be paired with
`arrivals = "join"` and an `arrivals.FUN` that returns those values,
`function(dat, at, new_ids, network) get_attr(dat, "hh_id")[new_ids]`,
so that the new nodes are connected to the members of their groups.

### Moving Nodes Between Groups

The clique edges are built from the grouping attribute once, when the
simulation starts. Setting the attribute with
[`set_attr()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
afterward does not rewire the layer: the old edges remain and the layer
no longer matches the groups.
[`move_to_group()`](https://epimodel.github.io/EpiModel/reference/move_to_group.md)
sets the attribute and the edges together, in either storage mode. Each
moved node loses its edges to the members of its old group and is
connected to every active member of its new group. In a model where
`age` is in years, a module in which young adults leave home to start
households of their own is:

``` r

leave_home <- function(dat, at) {
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  elig <- which(active == 1 & age >= 18 & age < 30)
  movers <- elig[runif(length(elig)) < 0.01]
  if (length(movers) > 0) {
    new_ids <- max(hh_id, na.rm = TRUE) + seq_along(movers)
    dat <- move_to_group(dat, ids = movers, group = new_ids)
  }
  return(dat)
}
```

The `ids` argument takes positional IDs. The `group` argument takes one
group ID per node, or a single value for all of them: an ID in use joins
that group, an unused ID starts a new group, and `NA` takes the node out
of its group. In a model with more than one clique layer, the `network`
argument names the layer to change.

The cumulative edgelist (see below) records the initial clique edges
with a `start` of 0. Edges added by arrivals and moves, and edges ended
by moves and departures, are recorded when the network resimulation
module next updates the cumulative edgelist, which it does for every
layer at every time step. A module that reads the cumulative edgelist
later in the same time step as a move should call
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
after the move.

## Observed Network Layers

[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
turns partial, usually egocentric, network data into a generative model
that
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
can simulate from. When the whole network has been observed (every node
and every contact, and for a dynamic network every time step), there is
nothing to estimate: the observed object is what
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
would otherwise have to generate.
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md)
wraps such a network, either a static `network` object or a
`networkDynamic` object with edge spells, as a layer that
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
accepts anywhere in its list of layers. Sensor, proximity-logger,
contact-tracing, and animal-tracking datasets are the usual sources.

Like a clique layer, an observed layer is an addition to
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
not a substitute for it. It is the wrong tool for a sample from which
one wants to generalize to a population, and it is not a way to model
how the observed ties arise; a network whose ties are to be reproduced
from their predictors is estimated with
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
whatever their turnover. The [Epidemics over Observed
Networks](https://epimodel.github.io/sismid/11_advanced/mod11-ObservedNets.html)
chapter of the NME course materials works through a complete example.

### Building an Observed Layer

``` r

library(networkDynamicData)
data(concurrencyComparisonNets)

obs <- netcensus(base)
obs
```

[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md)
carries the vertex attributes of the observed network into the
simulation as nodal attributes. Temporally extended vertex attributes
(those ending in `.active`), such as the `status.active` attribute
stored on `base`, are dropped with a message, since
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
keeps its own disease status. Printing a dynamic census shows the
observation window and the number of edges active per step across it;
printing a static census shows its edge count, mean degree, and number
of isolates. As with a clique layer, there is no model to diagnose, so
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) does
not accept a `netcensus` object.

### Simulation Time and the Observation Window

A static census has the same edges at every time step. For a dynamic
census, EpiModel time step `at` reads the edges active at time `at` of
the observed object, so the simulation clock is the observation clock.
The `window` argument of
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md)
sets the observation window. It defaults to the `net.obs.period` network
attribute when the object has one, and otherwise to the range of the
finite edge spell times; for `base`, the window runs from time 2 to time
102.

By the `networkDynamic` convention, edges active at the last observed
time stay active indefinitely, so a simulation that runs past the window
sees a frozen edge set. Keep `nsteps` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
within the window;
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
warns when `nsteps` reaches its end:

``` r

param <- param.net(inf.prob = 0.5, act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       resimulate.network = FALSE)
sim <- netsim(obs, param, init, control)
```

An observed layer is skipped by the network resimulation whatever the
value of `resimulate.network`. Setting it to `FALSE` is enough here
because the only layer is observed; the setting matters only when an
observed layer sits next to an estimated one.

### Observed Layers During Simulation

Without `tergmLite`, the observed `networkDynamic` object is used as
stored, and
[`get_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_edgelist.md)
reads the edges active at the current time step from it. Under
`tergmLite`, the edgelist of the layer is replaced at every time step
with the edges active in the observed object. Either way, the accessors
read the layer as any other, and the built-in infection modules transmit
across its edges. Next to other layers, the per-layer `inf.prob` and
`act.rate` are set through
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md):

``` r

# nw_obs: an observed static network on the same nodes as the estimated layer
sim <- netsim(list(netcensus(nw_obs), est),
              param.net(inf.prob = multilayer(0.3, 0.1), act.rate = 1),
              init, control)
```

An observed layer differs from an estimated one in a few other respects:

- **Fixed node set:** a census has a fixed node set, so
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  refuses vital dynamics in a model with an observed layer, and
  [`arrive_nodes()`](https://epimodel.github.io/EpiModel/reference/arrive_nodes.md)
  and
  [`depart_nodes()`](https://epimodel.github.io/EpiModel/reference/depart_nodes.md)
  stop if a custom module tries to add or remove nodes. Vertex activity
  spells are not used: every node is present throughout, and a node
  absent from part of the observation simply has no contacts then.
- **Network statistics:** with the default
  `nwstats.formula = "formation"`, the statistic recorded for the layer
  is its edge count, summarized at every time step for a dynamic census.
  There are no target statistics, so `print(sim)` and
  `plot(sim, type = "formation")` show the recorded values without
  targets.
- **Durations:** duration tracking under `tergmLite`
  (`tergmLite.track.duration = TRUE`) is refused for a dynamic census,
  whose observed spells already carry the edge durations.
- **Edge attributes:** edge attributes such as contact duration or
  contact count are not read, so the layer is binary. A custom module
  can read them from the observed object, which the layer’s network
  parameter record holds for a dynamic census as
  `get_nwparam(dat, network = 1)$census.nw`.
- **Module order:** a module that reads the layer before the network
  resimulation module runs in a time step sees the previous step’s edges
  under `tergmLite`, but the current step’s edges in full mode, which
  reads the `networkDynamic` object directly. Modules added to
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  under new names run before the built-in modules unless `module.order`
  places them later.

## Cumulative Edgelist

The cumulative edgelist is a historical record of all edges in a
network, including the time steps when each edge started and stopped.
This allows querying both current and past partnerships—essential for
contact tracing, partnership duration analysis, and reachability
analysis.

### Lifecycle

The cumulative edgelist follows a four-step lifecycle. The same data
structure is produced once and read at different stages:

1.  **Enable** in
    [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md):
    set `cumulative.edgelist = TRUE` (in-memory tracking during the run)
    and, if you want to use the result after
    [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
    returns, also set `save.cumulative.edgelist = TRUE`.
2.  **Track** during the run: the built-in network-resimulation module
    updates the cumulative edgelist once per network at every time step,
    using `control$truncate.el.cuml` to decide how much dissolved-edge
    history to retain. Custom modules that mutate the network outside
    the TERGM machinery should also call
    [`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
    after the change.
3.  **Read during the run** (inside a module) with
    [`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
    /
    [`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
    or with one of the derived helpers below
    ([`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
    [`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md)).
4.  **Read after the run** by accessing `sim$cumulative.edgelist[[s]]`
    directly on the returned `netsim` object.

The control settings, accessors, and helpers all share the same column
convention: `head`/`tail` (or `index`/`partner`) by unique ID,
`start`/`stop` time steps (inclusive, with `NA` `stop` for active
edges), and `network` index for multi-layer models.

### Enabling Cumulative Edgelists

Cumulative edgelist tracking must be explicitly enabled in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md):

``` r

control <- control.net(
  type = "SI",
  nsims = 1,
  nsteps = 100,
  cumulative.edgelist = TRUE,       # Enable in-memory tracking during the run
  truncate.el.cuml = 0,             # Drop dissolved edges immediately (default)
  save.cumulative.edgelist = TRUE,  # Attach the result to the returned netsim
  verbose = FALSE
)
```

Without `cumulative.edgelist = TRUE`, calls to
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
silently do nothing, and calls to
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
raise an error.

The `truncate.el.cuml` control sets the default truncation passed to
every automatic
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
call. The default (`0`) keeps only currently active edges, which is
enough for tracking active-edge start times while keeping memory low.
Use `Inf` to keep full history, or a positive integer to retain
dissolved edges for that many steps after they ended.

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

#### During Simulation: Specific Network

Inside a custom module, where `dat` is the live `netsim_dat` object:

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
12, inclusive. When no edges have been recorded yet, the function still
returns a `tibble` with these four columns and zero rows. If
`cumulative.edgelist = FALSE`, the function raises an error rather than
silently returning empty.

#### During Simulation: Multiple Networks

``` r

el_cumls <- get_cumulative_edgelists_df(dat, networks = NULL)
```

The `networks` argument accepts a vector of network indices or `NULL`
(all networks). The returned `data.frame` adds a `network` column
identifying which network layer each edge belongs to.

#### After Simulation

When `save.cumulative.edgelist = TRUE`, the cumulative edgelist is
attached to the `netsim` return object as `sim$cumulative.edgelist`, a
list with one element per simulation:

``` r

el_cuml <- sim$cumulative.edgelist[[1]]
```

Each element is the `data.frame` produced by
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md)
(`head`, `tail`, `start`, `stop`, `network`). Do not call the accessor
functions on a processed `netsim` object: the live `dat$run` state they
depend on is dropped during
[`process_out.net()`](https://epimodel.github.io/EpiModel/reference/process_out.net.md).
Note also that
[`truncate_sim()`](https://epimodel.github.io/EpiModel/reference/truncate_sim.md)
left-truncates the epidemiological time series but does not currently
trim `sim$cumulative.edgelist`; filter it manually by `start` / `stop`
if you need a matching time window.

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
`dat`—and are typically used for post-simulation analysis. Most often,
the input comes from the saved `sim$cumulative.edgelist[[s]]` slot of a
`netsim` run; see the *Converting External `networkDynamic` Objects* and
*Deduplicating Across Sources* subsections below for two ways to
construct input from other starting points.

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

### Converting External `networkDynamic` Objects

If you have a `networkDynamic` object that was produced outside of a
`cumulative.edgelist = TRUE` run (for example, from a non-`tergmLite`
simulation, from a saved RDS file, or built manually),
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md)
converts it into the same `data.frame` shape the reachability functions
expect:

``` r

nd <- get_network(sim, sim = 1, network = 1)   # full-mode netsim only
el_cuml <- as_cumulative_edgelist(nd)
fwd <- get_forward_reachable(el_cuml, from_step = 1, to_step = 100)
```

This is the recommended entry point for post-hoc reachability analysis
on simulations that did not enable cumulative-edgelist tracking up
front, provided the run kept the full `networkDynamic` history (i.e.,
not `tergmLite`).

### Deduplicating Across Sources

The reachability functions assume non-overlapping spells per
`(head, tail)` pair. When concatenating cumulative edgelists from
multiple sources (e.g., separate simulation segments, or several
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md)
conversions stitched together), use
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md)
to merge overlapping spells into a single row:

``` r

el_cuml <- dedup_cumulative_edgelist(rbind(el_run1, el_run2))
```

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
