# Working with Network Objects in EpiModel

## Introduction

This vignette discusses mechanisms usable inside `EpiModel` network
models with custom modules. More information about these in the
[Extending
EpiModel](https://epimodel.github.io/sismid/9_extending/mod9-Intro.html)
section of the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) course materials.

Inside the simulation, the networks themselves are stored under
`dat$run$nw`. Ultimately this vignette will describe multiple aspects of
working with network objects.

## Cumulative Edgelist

The cumulative edgelist refers to the historical list of edges in a
network with the time step they start and stopped. Such a list allows to
query current relationships (contacts, partnerships, etc.) as well as
past ones.

### Using the Cumulative Edgelist

The creation and update of the cumulative edgelist is done through the
[`EpiModel::update_cumulative_edgelist`](http://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
function.

``` r
dat <- update_cumulative_edgelist(dat, network, truncate = Inf)
```

This function takes 3 arguments:

1.  `dat`: the *Main List Object*.
2.  `network`: the number of the network for which to create the
    cumulative edgelist (for multi-layer network simulations).
3.  `truncate`: a number of time steps after which an inactive edge
    should be removed from the cumulative edgelist (this saves computer
    memory for large-scale simulations).

The function returns a modified version of `dat` that needs to be
assigned back.

The following snippet will update the cumulative edgelist for all the
networks used by a model and remove the edges that ended more than 100
steps ago.

``` r
for (n_network in seq_along(dat$run$nw)) {
  dat <- update_cumulative_edgelist(dat, n_network, truncate = 100)
}
```

In a complete model, this code would need to be run at the end of the
initialization module and at each time-step after the resimulation of
the networks.

### Accessing the Cumulative Edgelist

Cumulative edge-list refers to nodes with their Unique Ids. See
[`help("unique_id-tools", package = "EpiModel")`](http://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
for more information.

#### For a Specific Network

Accessing the cumulative edge-list of a given network is done using the
[`EpiModel::get_cumulative_edgelist`](http://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
function.

``` r
el_cuml <- get_cumulative_edgelist(dat, network)
```

The returned `el_cuml` object is a
[`tibble`](https://tibble.tidyverse.org/) with four columns:

1.  `head`: the `unique_id` first node of the edge.
2.  `tail`: the `unique_id` second node of the edge.
3.  `start`: the time-step where the edge was created.
4.  `stop`: the last time-step the edge was active.

[`EpiModel::get_cumulative_edgelist`](http://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
will **always** return a [`tibble`](https://tibble.tidyverse.org/) with
this 4 columns, even if the cumulative edgelist has not been calculated
for this particular network or if no edges are present. In these cases,
the [`tibble`](https://tibble.tidyverse.org/) will have no rows but keep
the correct column structure.

The `stop` column will **always** contain `NA` if an edge is currently
active.

Once an edge is not present anymore, the `stop` column for this edge
will contains the last step the edge was active. This means that an edge
with a `stop` value existed from `start` to `stop` both inclusive. This
makes it coherent with how `R` treats the indexes in a vector for
instance (from 1 to `length(vector)` inclusive).

#### For Multiple Networks

We often want to get the cumulative edgelist over several networks as
one.
[`EpiModel::get_cumulative_edgelists_df`](http://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md)
function provide such functionality.

``` r
el_cumls <- get_cumulative_edgelists_df(dat, networks = NULL)
```

The `networks` argument can be a vector of network position or `NULL`.
In this latter case, all networks will be selected.

The output of this function is similar to
[`EpiModel::get_cumulative_edgelist`](http://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
with the addition of a `network` column, indicating for each edge the
networks it exists on.

### Contact Tracing

A typical use of the cumulative edgelist is the trace the contacts of a
node over given number of steps. The
[`EpiModel::get_partners`](http://epimodel.github.io/EpiModel/reference/get_partners.md)
function simplifies this process:

``` r
partner_list <- get_partners(
    dat,
    index_posit_ids,
    networks = NULL,
    truncate = Inf,
    only.active.nodes = FALSE
)
```

Here we call “indexes” the nodes whose partners (contacts) we want to
extract. The arguments are:

1.  `dat`: as in `get_cumulative_edgelists_df`.
2.  `index_posit_ids`: a list of positional Ids for the indexes of
    interest.
3.  `networks`: as in `get_cumulative_edgelists_df`.
4.  `truncate`: similar to the `truncate` argument to
    `update_cumulative_edgelist` this argument filter out partnerships
    over this age.
5.  `only.active.nodes`: if set to `TRUE`, partnership with inactive
    nodes are removed.

The output is similar to `get_cumulative_edgelists_df` but the first two
columns are called `index` and `partner` and contains the Unique Ids of
the indexes given in argument in the first column and there partners in
the second one.

Note that the we refer to the indexes of interest with their Positional
Ids but the `index` and `partners` columns contains Unique Ids as they
can refer to nodes no longer in the network.
