# Functions to Access and Edit the Main netsim_dat Object in Network Models

These `get_`, `set_`, `append_`, and `add_` functions allow a safe and
efficient way to retrieve and mutate the main `netsim_dat` class object
of network models (typical variable name `dat`). They are intended for
use *inside module functions* that run during a
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
simulation, not for editing the user-facing `param.net`, `init.net`, or
`control.net` inputs before a run or the output object returned by
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).
See the **Intended Usage** section below for details and alternatives.

The accessors are documented on three pages, by the part of `dat` they
act on:

- **Nodal attributes**, with one value per node
  ([`get_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md)):
  [`get_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`get_attr_list()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`set_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`add_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`append_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`append_core_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  [`remove_node_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md),
  and
  [`get_core_attributes()`](https://epimodel.github.io/EpiModel/reference/get_attr.md).

- **Epidemic outputs**, with one value per time step
  ([`get_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md)):
  [`get_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md),
  [`get_epi_list()`](https://epimodel.github.io/EpiModel/reference/get_epi.md),
  [`set_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md),
  and
  [`add_epi()`](https://epimodel.github.io/EpiModel/reference/get_epi.md).

- **Parameters, control settings, and initial conditions**
  ([`get_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md)):
  [`get_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md),
  [`get_param_list()`](https://epimodel.github.io/EpiModel/reference/get_param.md),
  [`set_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md),
  [`add_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md),
  the same four functions for `control` and `init`, and
  [`get_network_control()`](https://epimodel.github.io/EpiModel/reference/get_param.md).

Nodes are referred to by their positional identifiers; see
[`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
and
[`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
to convert them to and from their unique identifiers.

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

## `set_` and `append_` vs `add_`

The `set_` and `append_` functions edit a pre-existing element or create
a new one if it does not exist already by calling the `add_` functions
internally.

## See also

[`update_params()`](https://epimodel.github.io/EpiModel/reference/update_params.md)
for editing a `param.net` object outside a simulation;
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
[`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md),
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
for input constructors;
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
for running a simulation.
