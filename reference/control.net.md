# Control Settings for Stochastic Network Models

Sets the controls for stochastic network models simulated with
[`netsim`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Usage

``` r
control.net(
  type,
  nsteps,
  start = 1,
  nsims = 1,
  ncores = 1,
  resimulate.network = FALSE,
  tergmLite = FALSE,
  cumulative.edgelist = FALSE,
  truncate.el.cuml = 0,
  attr.rules,
  epi.by,
  initialize.FUN = initialize.net,
  resim_nets.FUN = resim_nets,
  summary_nets.FUN = summary_nets,
  infection.FUN = NULL,
  recovery.FUN = NULL,
  departures.FUN = NULL,
  arrivals.FUN = NULL,
  nwupdate.FUN = nwupdate.net,
  prevalence.FUN = prevalence.net,
  verbose.FUN = verbose.net,
  module.order = NULL,
  save.nwstats = TRUE,
  nwstats.formula = "formation",
  save.transmat = TRUE,
  save.network,
  save.run = FALSE,
  save.cumulative.edgelist = FALSE,
  save.other,
  verbose = TRUE,
  verbose.int = 1,
  skip.check = FALSE,
  raw.output = FALSE,
  tergmLite.track.duration = FALSE,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 2e+05),
  set.control.tergm = control.simulate.formula.tergm(MCMC.maxchanges = Inf),
  save.diss.stats = TRUE,
  dat.updates = NULL,
  ...
)
```

## Arguments

- type:

  Disease type to be modeled, with the choice of `"SI"` for
  Susceptible-Infected diseases, `"SIR"` for
  Susceptible-Infected-Recovered diseases, and `"SIS"` for
  Susceptible-Infected-Susceptible diseases.

- nsteps:

  Number of time steps to simulate the model over. This must be a
  positive integer that is equal to the final step of a simulation. If a
  simulation is restarted with `start` argument, this number must be at
  least one greater than that argument's value.

- start:

  For models with network resimulation, time point to start up the
  simulation. For restarted simulations, this must be one greater than
  the final time step in the prior simulation and must be less than the
  value in `nsteps`.

- nsims:

  The total number of disease simulations.

- ncores:

  Number of processor cores to run multiple simulations on, using the
  `foreach` and `doParallel` implementations.

- resimulate.network:

  If `TRUE`, resimulate the network at each time step. This is required
  when the epidemic or demographic processes impact the network
  structure (e.g., vital dynamics).

- tergmLite:

  Logical indicating usage of either `tergm` (`tergmLite = FALSE`), or
  `tergmLite` (`tergmLite = TRUE`). Default of `FALSE`.

- cumulative.edgelist:

  If `TRUE`, calculates a cumulative edgelist within the network
  simulation module. This is used when tergmLite is used and the entire
  networkDynamic object is not used.

- truncate.el.cuml:

  Number of time steps of the cumulative edgelist to retain. See help
  for
  [`update_cumulative_edgelist`](http://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
  for options.

- attr.rules:

  A list containing the rules for setting the attributes of incoming
  nodes, with one list element per attribute to be set (see details
  below).

- epi.by:

  A character vector of length 1 containing a nodal attribute for which
  subgroup stratified prevalence summary statistics are calculated. This
  nodal attribute must be contained in the network model formation
  formula, otherwise it is ignored.

- initialize.FUN:

  Module to initialize the model at time 1, with the default function of
  [`initialize.net`](http://epimodel.github.io/EpiModel/reference/initialize.net.md).

- resim_nets.FUN:

  Module to resimulate the network at each time step, with the default
  function of
  [`resim_nets`](http://epimodel.github.io/EpiModel/reference/resim_nets.md).

- summary_nets.FUN:

  Module to extract summary statistics of the network at each time step,
  with the default function of
  [`summary_nets`](http://epimodel.github.io/EpiModel/reference/summary_nets.md).

- infection.FUN:

  Module to simulate disease infection, with the default function of
  [`infection.net`](http://epimodel.github.io/EpiModel/reference/infection.net.md).

- recovery.FUN:

  Module to simulate disease recovery, with the default function of
  [`recovery.net`](http://epimodel.github.io/EpiModel/reference/recovery.net.md).

- departures.FUN:

  Module to simulate departure or exit, with the default function of
  [`departures.net`](http://epimodel.github.io/EpiModel/reference/departures.net.md).

- arrivals.FUN:

  Module to simulate arrivals or entries, with the default function of
  [`arrivals.net`](http://epimodel.github.io/EpiModel/reference/arrivals.net.md).

- nwupdate.FUN:

  Module to handle updating of network structure and nodal attributes
  due to exogenous epidemic model processes, with the default function
  of
  [`nwupdate.net`](http://epimodel.github.io/EpiModel/reference/nwupdate.net.md).

- prevalence.FUN:

  Module to calculate disease prevalence at each time step, with the
  default function of
  [`prevalence.net`](http://epimodel.github.io/EpiModel/reference/prevalence.net.md).

- verbose.FUN:

  Module to print simulation progress to screen, with the default
  function of
  [`verbose.net`](http://epimodel.github.io/EpiModel/reference/verbose.net.md).

- module.order:

  A character vector of module names that lists modules in the order in
  which they should be evaluated within each time step. If `NULL`, the
  modules will be evaluated as follows: first any new modules supplied
  through `...` in the order in which they are listed, then the built-in
  modules in the order in which they are listed as arguments above.
  `initialize.FUN` will always be run first and `verbose.FUN` will
  always be run last.

- save.nwstats:

  If `TRUE`, save network statistics in a data frame. The statistics to
  be saved are specified in the `nwstats.formula` argument.

- nwstats.formula:

  A right-hand sided ERGM formula that includes network statistics of
  interest, with the default to the formation formula terms. Supports
  [`multilayer`](http://epimodel.github.io/EpiModel/reference/multilayer.md)
  specification.

- save.transmat:

  If `TRUE`, complete transmission matrix is saved at simulation end.

- save.network:

  If `TRUE`, networkDynamic or networkLite object is saved at simulation
  end.

- save.run:

  If `TRUE`, the `run` sublist of `dat` is saved, allowing a simulation
  to restart from this output.

- save.cumulative.edgelist:

  If `TRUE`, the `cumulative.edgelist` is saved at simulation end.

- save.other:

  A character vector of elements on the `netsim_dat` main data list to
  save out after each simulation. One example for base models is the
  attribute list, `"attr"`, at the final time step.

- verbose:

  If `TRUE`, print model progress to the console.

- verbose.int:

  Time step interval for printing progress to console, where `0` prints
  completion status of entire simulation and positive integer `x` prints
  progress after every `x` time steps. The default is to print progress
  after each time step.

- skip.check:

  If `TRUE`, skips the default error checking for the structure and
  consistency of the parameter values, initial conditions, and control
  settings before running base epidemic models. Setting this to `FALSE`
  is recommended when running models with new modules specified.

- raw.output:

  If `TRUE`, `netsim` will output a list of raw data (one per
  simulation) instead of a cleaned and formatted `netsim` object.

- tergmLite.track.duration:

  If `TRUE`, track duration information for models in `tergmLite`
  simulations. Supports
  [`multilayer`](http://epimodel.github.io/EpiModel/reference/multilayer.md)
  specification.

- set.control.ergm:

  Control arguments passed to
  [`ergm::simulate_formula.network`](https://rdrr.io/pkg/ergm/man/simulate.ergm.html).
  In `netsim`, this is only used when initializing the network with
  `edapprox = TRUE`. All other simulations in `netsim` use `tergm`.
  Supports
  [`multilayer`](http://epimodel.github.io/EpiModel/reference/multilayer.md)
  specification.

- set.control.tergm:

  Control arguments passed to
  [`tergm::simulate_formula.network`](https://rdrr.io/pkg/tergm/man/simulate.tergm.html).
  See the help file for
  [`netdx`](http://epimodel.github.io/EpiModel/reference/netdx.md) for
  details and examples on specifying this parameter. Supports
  [`multilayer`](http://epimodel.github.io/EpiModel/reference/multilayer.md)
  specification.

- save.diss.stats:

  If `TRUE`, `netsim` will compute and save duration and dissolution
  statistics for plotting and printing, provided `save.network` is
  `TRUE`, `tergmLite` is `FALSE`, and the dissolution model is
  homogeneous.

- dat.updates:

  Either `NULL`, a single function taking arguments `dat`, `at`, and
  `network`, or a list of functions of length one greater than the
  number of networks being simulated, with each function in the list
  taking arguments `dat` and `at`. Here `dat` is the main `netsim_dat`
  class object, `at` is the current timestep, and `network` is an index
  indicating the current position within the sequence of network
  (re)simulations on each time step. If a single function is passed, it
  will be called before the first network is simulated and after each
  network is simulated, with `network = 0L` before the first network is
  simulated and with `network = i` after the `i`th network is simulated.
  If a list of functions is passed, the first function will be called
  before the first network is simulated, and the `i + 1`th function will
  be called after the `i`th network is simulated. (Note that `at = 0L`
  is used for initial cross-sectional simulations in
  [`sim_nets_t1`](http://epimodel.github.io/EpiModel/reference/sim_nets_t1.md).)
  The function(s) should return the `netsim_dat` object with any updates
  needed to correctly represent the network states for calls to
  `simulate` and/or `summary`. This can be useful if nodal attributes
  appearing in one network model depend on nodal degrees in a different
  network.

- ...:

  Additional control settings passed to model.

## Value

An EpiModel object of class `control.net`.

## Details

`control.net` sets the required control settings for any network model
solved with the
[`netsim`](http://epimodel.github.io/EpiModel/reference/netsim.md)
function. Controls are required for both base model types and when
passing original process modules. For an overview of control settings
for base models, consult the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) course materials For all
base models, the `type` argument is a necessary parameter and it has no
default.

## The attr.rules Argument

The `attr.rules` parameter is used to specify the rules for how nodal
attribute values for incoming nodes should be set. These rules are only
necessary for models in which there are incoming nodes (i.e., arrivals).
There are three rules available for each attribute value:

- **`current`**: new nodes will be assigned this attribute in proportion
  to the distribution of that attribute among existing nodes at that
  current time step.

- **`t1`**: new nodes will be assigned this attribute in proportion to
  the distribution of that attribute among nodes at time 1 (that is, the
  proportions set in the original network for
  [`netest`](http://epimodel.github.io/EpiModel/reference/netest.md)).

- **`Value`**: all new nodes will be assigned this specific value, with
  no variation. For example, the rules list
  `attr.rules = list(race = "t1", sex = "current", status = "s")`
  specifies how the race, sex, and status attributes should be set for
  incoming nodes. By default, the rule is `"current"` for all attributes
  except status, in which case it is `"s"` (that is, all incoming nodes
  are susceptible).

## Checkpointing Simulations

`netsim` has a built-in checkpoint system to prevent losing computation
work if the function is interrupted (SIGINT, power loss, time limit
exceeded on a computation cluster). When enabled, each simulation will
be saved every `.checkpoint.steps` time steps. Then, if a checkpoint
enabled simulation is launched again with `netsim`, it will restart at
the last checkpoint available in the saved data.

To enable the checkpoint capabilities of `netsim`, two control arguments
have to be set: `.checkpoint.steps`, which is a positive number of time
steps to be run between each file save; and `.checkpoint.dir`, which is
the path to a directory to save the checkpointed data. If
`.checkpoint.dir` directory does not exist, `netsim` will attempt to
create it on the first checkpoint save. With these two controls defined,
one can simply re-run `netsim` with the same arguments to restart a set
of simulations that were interrupted.

Simulations are checkpointed individually: for example, if 3 simulations
are run on a single core, the first 2 are finished, then the
interruption occurs during the third, `netsim` will only restart the
third one from the last checkpoint.

A `.checkpoint.compress` argument can be set to overwrite the `compress`
argument in `saveRDS` used to save the checkpointed data. The current
default for `saveRDS` is `gunzip (gz)`, which provides fast compression
that usually works well on `netsim` objects.

By default, if `netsim` reaches the end of all simulations, the
checkpoint data directory and its content are removed before returning
the `netsim` object. The `.checkpoint.keep` argument can be set to
`TRUE` to prevent this removal to inspect the raw simulation objects.

## New Modules

Base network models use a set of module functions that specify how the
individual nodes in the network are subjected to infection, recovery,
demographics, and other processes. Core modules are those listed in the
`.FUN` arguments. For each module, there is a default function used in
the simulation. The default infection module, for example, is contained
in the
[`infection.net`](http://epimodel.github.io/EpiModel/reference/infection.net.md)
function.

For original models, one may substitute replacement module functions for
any of the default functions. New modules may be added to the workflow
at each time step by passing a module function via the `...` argument.
Consult the [Extending
EpiModel](https://epimodel.github.io/sismid/9_extending/mod9-Intro.html)
section of the Network Modeling for Epidemics course materials. One may
remove existing modules, such as `arrivals.FUN`, from the workflow by
setting the parameter value for that argument to `NULL`.

## End Horizon

`netsim` implements an "End Horizon" mechanism, where a set of modules
are removed from the simulation at a specific time step. This is enabled
through the `end.horizon` parameter to `control.net`.

This parameter must receive a `list` with fields `at`, the time step at
which the end horizon occurs, and `modules`, a character vector with the
names of the modules to remove. (e.g \`list(at = 208, modules =
c("arrivals.FUN", "infections.FUN")))

## See also

Use
[`param.net`](http://epimodel.github.io/EpiModel/reference/param.net.md)
to specify model parameters and
[`init.net`](http://epimodel.github.io/EpiModel/reference/init.net.md)
to specify the initial conditions. Run the parameterized model with
[`netsim`](http://epimodel.github.io/EpiModel/reference/netsim.md).
