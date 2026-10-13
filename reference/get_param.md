# Get and Set Parameters, Controls, and Initial Conditions in Network Models

Retrieve, create, and edit the parameters, control settings, and initial
conditions of the main `netsim_dat` object of network models, from
inside module functions. They start as the
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md),
and
[`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md)
inputs of
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md),
and the three share one set of accessors:

- `get_param()`, `get_control()`, and `get_init()` return one element.
  `get_param_list()`, `get_control_list()`, and `get_init_list()` return
  several elements as a named list.

- `set_param()`, `set_control()`, and `set_init()` replace the value of
  an element, and create the element if it does not exist.
  `add_param()`, `add_control()`, and `add_init()` create an element set
  to `NA`.

- `get_network_control()` returns the entry for one network layer of a
  control setting given as a
  [`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
  object.

Unlike
[`set_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md)
and
[`set_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md),
these `set_` functions take neither nodes nor a time step: the new value
is used from then on, until it is set again.

See
[net-accessor](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
for an overview of the accessors.

## Usage

``` r
get_param(dat, item, override.null.error = FALSE)

get_param_list(dat, item = NULL)

set_param(dat, item, value)

add_param(dat, item)

get_control(dat, item, override.null.error = FALSE)

get_control_list(dat, item = NULL)

get_network_control(dat, network, item, override.null.error = FALSE)

set_control(dat, item, value)

add_control(dat, item)

get_init(dat, item, override.null.error = FALSE)

get_init_list(dat, item = NULL)

set_init(dat, item, value)

add_init(dat, item)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  Name of the parameter, control setting, or initial condition. For the
  `get_*_list()` functions, a character vector of names, or `NULL`
  (default) for all elements.

- override.null.error:

  If `TRUE`, return `NULL` instead of an error when `item` does not
  exist (default = `FALSE`).

- value:

  New value of the element.

- network:

  Index of the network layer.

## Value

For `get_param()`, `get_control()`, and `get_init()`, the value of the
element; for `get_network_control()`, its entry for layer `network`; for
the `get_*_list()` functions, a named list of elements. The other
functions return the updated `netsim_dat` object.

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
[`get_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md)
and
[`get_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md)
for the other groups;
[`update_params()`](https://epimodel.github.io/EpiModel/reference/update_params.md)
for editing a `param.net` object outside a simulation.

## Examples

``` r
dat <- create_dat_object(control = list(nsteps = 150))

dat <- add_param(dat, "x")
get_param(dat, "x")
#> [1] NA
dat <- set_param(dat, "x", 0.4)
dat <- set_param(dat, "y", 0.8)
get_param_list(dat)
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_param(dat, "x")
#> [1] 0.4

# The control and init accessors work the same way
dat <- add_control(dat, "x")
dat <- set_control(dat, "x", 0.4)
get_control_list(dat, c("nsteps", "x"))
#> $nsteps
#> [1] 150
#> 
#> $x
#> [1] 0.4
#> 
get_control(dat, "nsteps")
#> [1] 150

dat <- add_init(dat, "x")
dat <- set_init(dat, "x", 0.4)
get_init_list(dat)
#> $x
#> [1] 0.4
#> 
get_init(dat, "x")
#> [1] 0.4

# A control setting with one entry per network layer
dat <- set_control(dat, "z", multilayer(0.1, 0.2))
get_network_control(dat, network = 2, "z")
#> [1] 0.2
```
