# Observed Network Layer for Network Epidemic Models

Wraps a fully observed contact network, either a static `network` or a
`networkDynamic` with edge spells, as a layer that
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
accepts anywhere in its layer list alongside the layers estimated with
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md).
The simulation reads the edges active at each time step from the
observed object; there is no formation model, no dissolution model, and
no resimulation. It is an addition to `netest` for the case in which the
whole network has been observed, not a substitute for estimation from a
sample (see Details).

## Usage

``` r
netcensus(nw, window = NULL)
```

## Arguments

- nw:

  An object of class `network` (a static census: the edges are active at
  every time step) or `networkDynamic` (a dynamic census: the edges
  active at time step `at` are those whose spells cover `at`). The
  network must be undirected and not bipartite. Vertex attributes on
  `nw` are carried into the simulation as nodal attributes; temporally
  extended vertex attributes (`*.active`) are dropped with a message,
  since `netsim` keeps its own disease status.

- window:

  Observation window of a dynamic census, as `c(start, end)`. Defaults
  to the `net.obs.period` network attribute when `nw` has one, and
  otherwise to the range of the finite edge spell times. Used to warn
  when a simulation runs past the end of the observations.

## Value

An object of class `netcensus`, a list with elements:

- **newnetwork:** the observed network, with temporally extended vertex
  attributes removed.

- **dynamic:** `TRUE` for a `networkDynamic` census, `FALSE` for a
  static one.

- **window:** the observation window of a dynamic census, or `NULL`.

- **edapprox:** `FALSE`.

- **summary:** a list with the node count, the vertex attribute names,
  and either the edge count, mean degree, and number of isolates
  (static) or the number of distinct edges ever observed and the number
  of edges active at each integer time in the window (dynamic).

## Details

[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
exists to turn partial, usually egocentric, network data into a
generative model that `netsim` can simulate from. That is the general
workflow, because samples are the norm and the fit is what lets a sample
stand in for a population. When the whole network has been observed, for
every node and every contact and, for a dynamic census, over every time
step, there is nothing to estimate: the observed object is what `netsim`
would otherwise have to generate. `netcensus` places it in the layer
list directly. Sensor, proximity-logger, contact-tracing, and
animal-tracking datasets are the usual sources. It is the wrong tool for
a sample from which one wants to generalize, and it is not a way to
model how the observed ties arise; a network whose ties are to be
reproduced from their predictors is estimated with `netest`, whatever
their turnover.

For a `networkDynamic` census, EpiModel time step `at` reads the edges
active at time `at` of the observed object, so the simulation clock is
the observation clock. Match `nsteps` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
to the observation window: by the `networkDynamic` convention, edges
active at the last observed time stay active indefinitely, so a
simulation that runs past the window sees a frozen edge set.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
warns when `nsteps` reaches the end of `window`. Vertex activity spells
are not used; every node is present throughout.

A census has a fixed node set, so vital dynamics are refused:
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
stops when the model has arrivals or departures and a `netcensus` layer,
and
[`arrive_nodes()`](https://epimodel.github.io/EpiModel/reference/arrive_nodes.md)
and
[`depart_nodes()`](https://epimodel.github.io/EpiModel/reference/depart_nodes.md)
stop if a custom module tries to add or remove nodes. Duration tracking
under `tergmLite` is also refused for a dynamic census, whose observed
spells already carry the durations.

The layer is skipped by the network resimulation and the edges
correction, and the built-in infection modules read it through
[`discord_edgelist()`](https://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
like any other layer, so per-layer `inf.prob` and `act.rate` may be
given with
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
in
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md).
Under `tergmLite`, the active edgelist is extracted from the observed
object at every time step; without `tergmLite`, the `networkDynamic`
object is used as stored. Network statistics are recorded through
`nwstats.formula` in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md),
with the default `"formation"` meaning `~edges`.
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
refuses a `netcensus` object, since there is no model to diagnose;
`print` summarizes the observed edges over the window.

## See also

[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
for a layer of cliques defined by a grouping attribute, which is the
other model-free layer.
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
runs the simulation, and
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
sets per-layer parameters and controls.

## Examples

``` r
if (FALSE) { # \dontrun{
# A dynamic census: an observed networkDynamic with edge spells
library(networkDynamicData)
data(concurrencyComparisonNets)
obs <- netcensus(base)
obs

param <- param.net(inf.prob = 0.5, act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       resimulate.network = FALSE, verbose = FALSE)
sim <- netsim(obs, param, init, control)
plot(sim)

# The same census under tergmLite, with the active edges read each step
control <- control.net(type = "SI", nsteps = 100, nsims = 5,
                       tergmLite = TRUE, verbose = FALSE)
sim <- netsim(obs, param, init, control)

# A static census next to an estimated layer
nw <- network_initialize(n = 100)
nw <- add.edges(nw, tail = 1:50, head = 51:100)
est <- netest(nw, ~edges, target.stats = 40,
              coef.diss = dissolution_coefs(~offset(edges), 10),
              verbose = FALSE)
sim <- netsim(list(netcensus(nw), est),
              param.net(inf.prob = multilayer(0.3, 0.1)),
              init.net(i.num = 5),
              control.net(type = "SI", nsteps = 20, nsims = 1,
                          tergmLite = TRUE, verbose = FALSE))
} # }
```
