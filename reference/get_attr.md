# Get and Set Nodal Attributes in Network Models

Retrieve, create, and edit the nodal attributes of the main `netsim_dat`
object of network models, from inside module functions. Each nodal
attribute is a vector with one value per node, indexed by the positional
identifiers of the nodes (`posit_ids`).

- `get_attr()` returns one attribute, for all nodes or for the nodes in
  `posit_ids`. `get_attr_list()` returns several attributes as a named
  list.

- `set_attr()` replaces the values of an attribute, for all nodes or for
  the nodes in `posit_ids`, and creates the attribute if it does not
  exist. `add_attr()` creates an attribute set to `NA` for all nodes.

- `append_core_attr()` adds `n.new` nodes and sets their core
  attributes. `append_attr()` then adds the values of each other
  attribute for the same nodes.

- `remove_node_attr()` removes the nodes in `posit_ids` from every
  attribute.

- `get_core_attributes()` lists the core attributes with their types.

See
[net-accessor](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
for an overview of the accessors.

## Usage

``` r
get_attr(dat, item, posit_ids = NULL, override.null.error = FALSE)

get_attr_list(dat, item = NULL)

set_attr(dat, item, value, posit_ids = NULL, override.length.check = FALSE)

add_attr(dat, item)

append_attr(dat, item, value, n.new)

append_core_attr(dat, at, n.new)

remove_node_attr(dat, posit_ids)

get_core_attributes()
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  Name of the nodal attribute. For `get_attr_list()`, a character vector
  of names, or `NULL` (default) for all attributes.

- posit_ids:

  Positional identifiers of the nodes. For `get_attr()` and
  `set_attr()`, the nodes to get or edit, with `NULL` (default) for all
  nodes; for `remove_node_attr()`, the nodes to remove.

- override.null.error:

  If `TRUE`, return `NULL` instead of an error when `item` does not
  exist (default = `FALSE`).

- value:

  New values of the attribute. For `set_attr()`, one value per node, or
  with `posit_ids`, one value per selected node or a single value for
  all of them. For `append_attr()`, one value per new node or a single
  value for all of them.

- override.length.check:

  If `TRUE`, `set_attr()` accepts a `value` whose length differs from
  the number of nodes (default = `FALSE`). Every attribute must again
  hold one value per node by the end of the time step.

- n.new:

  Number of new nodes: for `append_core_attr()`, the nodes to add; for
  `append_attr()`, the number of values to append at the end of `item`.

- at:

  Current time step, recorded as the entry time (`entrTime`) of the new
  nodes.

## Value

For `get_attr()`, the attribute vector, or its values for the nodes in
`posit_ids`; for `get_attr_list()`, a named list of attribute vectors;
for `get_core_attributes()`, a named list with a value of the correct
type for each core attribute. The other functions return the updated
`netsim_dat` object.

## Core Attributes

The `append_core_attr` function initializes the attributes necessary for
EpiModel to work (the four core attributes are: "active", "unique_id",
"entrTime", and "exitTime"). These attributes are used in the
initialization phase of the simulation, to create the nodes (see
[`initialize.net()`](https://epimodel.github.io/EpiModel/reference/initialize.net.md));
and also used when adding nodes during the simulation (see
[`arrivals.net()`](https://epimodel.github.io/EpiModel/reference/arrivals.net.md)).

## One Value Per Node

Every nodal attribute holds exactly one value per node, and
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
checks this at the end of each time step. When `set_attr` and
`append_attr` create an attribute, the nodes they do not write to are
filled with `NA` to preserve that: `set_attr` with `posit_ids` fills the
unselected nodes, and `append_attr` fills the nodes that predate the
ones being appended. This means `append_core_attr` must be called before
`append_attr` when adding new nodes, as it is what makes those nodes
known to the rest of the object.

## Intended Usage

These accessors operate on the live `netsim_dat` object that EpiModel
constructs internally when a simulation starts and passes through each
module at every time step. The expected calling context is **inside a
user-supplied module function** (e.g., a custom `infection.FUN`,
`recovery.FUN`, or arrivals/departures module) registered with
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).

They are *not* a general-purpose tool for editing the user-facing input
or output objects of
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md):

- Calling them on a
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
  [`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md),
  or
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  object before a run will appear to succeed (those objects are also
  list-like) but produces an object that is not a valid simulation
  input.

- Calling them on the object returned by
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  modifies only the saved input slot, not anything that would change the
  result of a re-run. Note also that
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  strips the `*.FUN` entries from its returned `control` slot, so
  feeding that slot back into a new simulation will fail.

For modifications outside a running simulation, use the following
instead:

- **Editing parameters before a run:**
  [`update_params()`](https://epimodel.github.io/EpiModel/reference/update_params.md)
  for a `param.net` object, or direct list assignment (e.g.,
  `p$inf.prob <- 0.5`).

- **Editing init or control before a run:** direct list assignment
  (e.g., `ctrl$nsteps <- 1000`), or rebuild with a fresh call to the
  [`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md)
  /
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  constructor.

- **Scheduled mid-simulation changes:** the `.param.updater.list`
  argument to
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
  and the `.control.updater.list` argument to
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).
  See the "model-parameters" vignette for the underlying updater module
  and the scenario API built on top of it.

## Mutability

The `set_`, `append_`, and `add_` functions DO NOT modify the
`netsim_dat` object in place. The result must be assigned back to `dat`
in order to be registered: `dat <- set_*(dat, item, value)`.

## See also

[net-accessor](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
for an overview of the accessors, and
[`get_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md)
and
[`get_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md)
for the other groups;
[`arrive_nodes()`](https://epimodel.github.io/EpiModel/reference/arrive_nodes.md)
and
[`depart_nodes()`](https://epimodel.github.io/EpiModel/reference/depart_nodes.md)
for adding and removing nodes from the networks;
[`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
and
[`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
for converting node identifiers.

## Examples

``` r
dat <- create_dat_object(control = list(nsteps = 150))
dat <- append_core_attr(dat, at = 1, n.new = 10)

dat <- add_attr(dat, "age")
dat <- set_attr(dat, "age", runif(10, 15, 65))
dat <- set_attr(dat, "status", rep(c("s", "i"), c(8, 2)))
get_attr_list(dat)
#> $active
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
#> $entrTime
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
#> $exitTime
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> $unique_id
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $age
#>  [1] 55.47115 35.53215 46.12376 19.40064 42.05494 45.12422 59.13545 24.87035
#>  [9] 44.30690 17.36665
#> 
#> $status
#>  [1] "s" "s" "s" "s" "s" "s" "s" "s" "i" "i"
#> 
get_attr_list(dat, c("age", "active"))
#> $age
#>  [1] 55.47115 35.53215 46.12376 19.40064 42.05494 45.12422 59.13545 24.87035
#>  [9] 44.30690 17.36665
#> 
#> $active
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
get_attr(dat, "status")
#>  [1] "s" "s" "s" "s" "s" "s" "s" "s" "i" "i"
get_attr(dat, "status", posit_ids = c(1, 9))
#> [1] "s" "i"

# Edit selected nodes
dat <- set_attr(dat, "status", "i", posit_ids = 3)
get_attr(dat, "status")
#>  [1] "s" "s" "i" "s" "s" "s" "s" "s" "i" "i"

# Add 5 nodes: core attributes first, then the others
dat <- append_core_attr(dat, at = 2, n.new = 5)
dat <- append_attr(dat, "status", "s", 5)
dat <- append_attr(dat, "age", 15, 5)
get_attr(dat, "age")
#>  [1] 55.47115 35.53215 46.12376 19.40064 42.05494 45.12422 59.13545 24.87035
#>  [9] 44.30690 17.36665 15.00000 15.00000 15.00000 15.00000 15.00000
get_attr(dat, "entrTime")
#>  [1] 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
```
