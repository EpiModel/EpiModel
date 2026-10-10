# Clique Layer for Network Epidemic Models

Builds a network layer in which every pair of nodes sharing a value of a
grouping attribute is connected, so that each group (a household,
classroom, hospital ward, or cabin) is a clique.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
accepts the layer anywhere in its layer list alongside the layers
estimated with
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md).
The edges of a clique layer never form or dissolve: there is no
formation model, no dissolution model, and no resimulation. It is an
addition to `netest` for structures that are groups by definition, not a
substitute for an ERGM (see Details).

## Usage

``` r
netclique(
  nw,
  group.attr,
  arrivals = c("isolate", "new", "join"),
  arrivals.FUN = NULL
)
```

## Arguments

- nw:

  An object of class `network` or `networkLite` holding the node set and
  any vertex attributes, as passed to
  [`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
  for the other layers of the same model. The network must be undirected
  and not bipartite. Any edges on `nw` are ignored.

- group.attr:

  Name of a vertex attribute on `nw` whose values partition the nodes
  into groups. Every pair of nodes sharing a non-missing value becomes
  an edge. Nodes with a missing (`NA`) value are in no group and are
  isolates on this layer.

- arrivals:

  Rule for placing nodes that arrive during a simulation with vital
  dynamics: `"isolate"` (the default; the new node has no edges on this
  layer and its group attribute is `NA`), `"new"` (each arrival starts a
  group of size one, with a fresh group id), or `"join"` (each arrival
  joins an existing group and is connected to all of its current
  members). See Details.

- arrivals.FUN:

  Optional function used with `arrivals = "join"` to choose the group
  each arrival joins. See Details.

## Value

An object of class `netclique`, a list with elements:

- **newnetwork:** a `networkLite` object holding the node set, the
  vertex attributes of `nw`, and the clique edges. A `networkLite`
  stores the edges as an edgelist, so the layer stays small for large
  populations;
  [`trim_netest()`](https://epimodel.github.io/EpiModel/reference/trim_netest.md)
  has nothing further to remove from it.

- **group.attr:** the grouping attribute name.

- **arrivals**, **arrivals.FUN:** the arrival rule.

- **target.stats**, **target.stats.names:** the edge count of the layer,
  named `"edges"`, so that the network statistics table printed by
  [`print.netsim()`](https://epimodel.github.io/EpiModel/reference/print.netsim.md)
  and plotted by
  [`plot.netsim()`](https://epimodel.github.io/EpiModel/reference/plot.netsim.md)
  shows the layer's own edge count as the target.

- **edapprox:** `FALSE`.

- **summary:** a list with the node count, number of groups, edge count,
  mean degree, number of isolates, and the group size distribution.

## Details

A `netest` object carries a formation model, a dissolution model, and a
starting network; `netsim` resimulates its layer each time step. That is
the right representation whenever tie existence is something to model: a
network with a degree distribution, mixing by attribute, and clustering
to reproduce is an ERGM's job whether its ties turn over quickly,
slowly, or not at all, and a long partnership duration in
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
keeps a `netest` layer close to fixed.

Some contact structures are groups by construction instead: every pair
of co-residents is a household contact for the whole simulation, every
pair of pupils in a classroom is a classroom contact, every pair of
cabin-mates shares a cabin. There is no tie-formation process to
estimate, and an ERGM can reproduce such cliques only through a
`nodematch` term on the group attribute targeted at its maximum, the
number of within-group pairs. A target at the maximum of a statistic has
no finite coefficient, so the fit runs to its iteration limit, the
coefficient it stops at is arbitrary, and a few groups are typically
left incomplete. Before this function, the alternative was to carry the
clique edgelist into `netsim` as a parameter and walk it in a custom
infection module. `netclique` is for these cases. It builds a layer
object that fills a slot in the layer list, is skipped by the network
resimulation and the edges correction, and is read by the built-in
infection modules through
[`discord_edgelist()`](https://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
like any other layer.

A network whose edges are given rather than generated, such as a fully
observed contact network, is a different case again: it may carry its
own edge dynamics and it skips estimation for a different reason. That
case is handled by
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md).

The group attribute must be an integer, numeric, or character vector; a
factor is refused because arrivals under `"new"` need to create values
the factor does not have. It should also be set on the network passed to
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
for the other layers, since `netsim` reads nodal attributes from the
first layer in its list; when it is not, `netsim` copies it from the
clique layer.

## Arrivals

When a model has vital dynamics, each arriving node has to be placed on
the clique layer. `netsim` applies the rule chosen with `arrivals`
inside
[`arrive_nodes()`](https://epimodel.github.io/EpiModel/reference/arrive_nodes.md),
once per clique layer, after the arrivals module has created the node
and set its other attributes:

- `"isolate"`: the node gets no edges and its group attribute is `NA`. A
  custom module may place it later.

- `"new"`: the node gets a fresh group id and no edges. Fresh ids
  continue from the largest id used so far in the simulation, so an id
  is never given to a second group, even after every member of the first
  has departed; for character ids they are new strings of the form
  `arrival_<i>`. Groups then grow only through further `"join"`
  arrivals, so this rule suits models in which arrivals are
  single-person households.

- `"join"`: the node joins an existing group and is connected to every
  active member of that group, including other nodes joining the same
  group in the same time step. By default the group is drawn with
  probability proportional to its current size, which is the
  distribution of a randomly chosen existing member. `arrivals.FUN`
  replaces that draw: it is called as
  `arrivals.FUN(dat, at, new_ids, network)` and must return one group id
  per element of `new_ids` (an `NA` leaves that node as an isolate).
  Because the function receives the full `dat` object it can read any
  nodal attribute; the example below places each newborn in a household
  that already has a young child, and a function that reads the group
  attribute back for `new_ids` lets a custom arrivals module set the
  group itself. A function that draws from existing groups should leave
  `new_ids` and inactive nodes out of its pool: when an arrivals module
  gives new nodes a default value of the grouping attribute, such as
  `0`, a pool taken from all nodes contains that value, and every
  arrival that draws it joins the other arrivals in a spurious group.

The layer keeps `arrivals.FUN` with its enclosing environment, and so
does every `netsim` object simulated from it. Define the function at the
top level of a script or in a package, not inside another function whose
local objects (fitted networks, population data) would then be saved
with each simulation.

The rule decides the group attribute of every arriving node. A value
given by `attr.rules` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
or by a custom arrivals module is replaced, so that the layer always
matches the groups, and `netsim` warns when `attr.rules` names the group
attribute. The exception is `"join"` with an `arrivals.FUN` that reads
that value back, as described above.

## Transmission over a clique layer

The built-in infection modules treat every layer alike, so a clique
layer transmits with the model's `inf.prob` and `act.rate` unless those
are given per layer as
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
objects in
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md):
`param.net(inf.prob = multilayer(0.45, 0.10), act.rate = multilayer(1, 2))`
with the household layer first and the community layer second. The
transmission matrix records the layer of each transmission in its
`network` column.

Nodes departing the population are removed from a clique layer with
their edges, as on every other layer. Network statistics for the layer
are recorded through `nwstats.formula` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md),
with the default `"formation"` meaning `~edges` for a clique layer, and
the cumulative edgelist records its edges with a start time of 0.
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
refuses a `netclique` object, since there is no model to diagnose; the
group size distribution and mean degree are shown by `print`.

## See also

[`sample_groups()`](https://epimodel.github.io/EpiModel/reference/sample_groups.md)
builds a population from a table of group types, and
[`assign_groups()`](https://epimodel.github.io/EpiModel/reference/assign_groups.md)
assigns group ids to an existing population with composition
constraints.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
runs the simulation, and
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
sets per-layer parameters and controls.

## Examples

``` r
if (FALSE) { # \dontrun{
# Households as cliques under a community TERGM layer
hh <- sample_groups(500, c("adult" = 0.3, "adult adult" = 0.4,
                           "adult adult child" = 0.2,
                           "adult adult child child" = 0.1),
                    attr.name = "age")
nw <- network_initialize(n = 500)
nw <- set_vertex_attribute(nw, "age", hh$age)
nw <- set_vertex_attribute(nw, "hh_id", hh$group)

est_hh <- netclique(nw, group.attr = "hh_id", arrivals = "join")
est_hh
print(est_hh, by = "age")

est_com <- netest(nw, formation = ~edges + nodematch("age"),
                  target.stats = c(250, 175),
                  coef.diss = dissolution_coefs(~offset(edges), 20),
                  verbose = FALSE)

param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 50, nsims = 1,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       verbose = FALSE)
sim <- netsim(list(est_hh, est_com), param, init, control)

# Share of transmissions on the household layer
tm <- get_transmat(sim)
mean(tm$network == 1)

# Newborns join a household that already has a child
est_hh2 <- netclique(nw, group.attr = "hh_id", arrivals = "join",
  arrivals.FUN = function(dat, at, new_ids, network) {
    active <- get_attr(dat, "active")
    hh <- get_attr(dat, "hh_id")
    age <- get_attr(dat, "age")
    existing <- active == 1
    existing[new_ids] <- FALSE
    pool <- hh[which(existing & age == "child" & !is.na(hh))]
    pool[sample.int(length(pool), length(new_ids), replace = TRUE)]
  })
} # }
```
