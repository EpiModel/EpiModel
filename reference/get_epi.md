# Get and Set Epidemic Outputs in Network Models

Retrieve, create, and edit the epidemic outputs of the main `netsim_dat`
object of network models, from inside module functions. Each epidemic
output is a vector with one value per time step, which becomes a column
of the data frame returned by
[`as.data.frame.netsim()`](https://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md)
after the simulation.

- `get_epi()` returns one output, for all time steps or for the time
  steps in `at`. `get_epi_list()` returns several outputs as a named
  list.

- `set_epi()` records `value` at one time step, the current one by
  default, and creates the output if it does not exist. `add_epi()`
  creates an output set to `NA` at every time step.

Unlike the other `set_` functions, `set_epi()` takes a time step `at`,
which comes last and must be passed by name.

See
[net-accessor](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
for an overview of the accessors.

## Usage

``` r
get_epi(dat, item, at = NULL, override.null.error = FALSE)

get_epi_list(dat, item = NULL)

set_epi(dat, item, value, ..., at = NULL)

add_epi(dat, item)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  Name of the epidemic output. For `get_epi_list()`, a character vector
  of names, or `NULL` (default) for all outputs.

- at:

  For `get_epi()`, the time steps at which to return the output, with
  `NULL` (default) for all time steps. For `set_epi()`, the time step at
  which to record `value`, which defaults to the current time step
  ([`get_current_timestep()`](https://epimodel.github.io/EpiModel/reference/get_current_timestep.md))
  when `NULL`.

- override.null.error:

  If `TRUE`, return `NULL` instead of an error when `item` does not
  exist (default = `FALSE`).

- value:

  Value of the output at time step `at`.

- ...:

  Must be empty. It makes `at` a named-only argument of `set_epi()`, so
  that calls written for the `set_epi(dat, item, at, value)` signature
  of EpiModel 2.6 and earlier error instead of silently swapping `at`
  and `value`.

## Value

For `get_epi()`, the output vector, or its values at the time steps in
`at`; for `get_epi_list()`, a named list of output vectors. The other
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
[`get_param()`](https://epimodel.github.io/EpiModel/reference/get_param.md)
for the other groups;
[`epi_trackers()`](https://epimodel.github.io/EpiModel/reference/epi_trackers.md)
for outputs computed by tracker functions at each time step;
[`mutate_epi()`](https://epimodel.github.io/EpiModel/reference/mutate_epi.md)
for outputs derived after the simulation.

## Examples

``` r
dat <- create_dat_object(control = list(nsteps = 5))
dat <- set_current_timestep(dat, 1)

dat <- add_epi(dat, "i.num")
dat <- set_epi(dat, "i.num", 10)
dat <- set_epi(dat, "s.num", 90)
get_epi_list(dat)
#> $i.num
#> [1] 10 NA NA NA NA
#> 
#> $s.num
#> [1] 90 NA NA NA NA
#> 

# Record at a time step other than the current one
dat <- set_epi(dat, "i.num", 12, at = 2)
get_epi_list(dat, c("i.num", "s.num"))
#> $i.num
#> [1] 10 12 NA NA NA
#> 
#> $s.num
#> [1] 90 NA NA NA NA
#> 
get_epi(dat, "i.num")
#> [1] 10 12 NA NA NA
get_epi(dat, "i.num", at = c(1, 2))
#> [1] 10 12
```
