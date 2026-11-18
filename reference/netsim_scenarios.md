# Run EpiModel Network Simulations With Scenarios

This function will run `n_rep` replications of each scenarios in the
`scenarios_list`. It runs them as multiple batches of up to `n_cores`
simulations at a time. The simfiles are then stored in the `output_dir`
folder and are named using the following pattern:
"sim\_\_name_of_scenario\_\_2.rds". Where the last number is the batch
number for this particular scenario. Each scenario is therefore run over
`ceiling(n_rep / n_cores)` batches.

## Usage

``` r
netsim_scenarios(
  path_to_x,
  param,
  init,
  control,
  scenarios_list,
  n_rep,
  n_cores,
  output_dir,
  ...
)
```

## Arguments

- path_to_x:

  Path to a fitted network model or a restart point saved with `saveRDS`
  (See the `x` argument to the
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  function), or to a restart pool directory (see
  [`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md)).

- param:

  Model parameters, as an object of class
  [`param.net`](https://epimodel.github.io/EpiModel/reference/param.net.md).
  Includes transmission probability (`inf.prob`), act rate (`act.rate`),
  recovery rate (`rec.rate`), and demographic rates for models with
  vital dynamics. Custom parameters for extended models may also be
  passed through `param.net`.

- init:

  Initial conditions, as an object of class
  [`init.net`](https://epimodel.github.io/EpiModel/reference/init.net.md).
  Specifies the initial number of infected nodes (`i.num`) and, for SIR
  models, recovered nodes (`r.num`). For two-group models, the
  corresponding `.g2` parameters are also required.

- control:

  Control settings, as an object of class
  [`control.net`](https://epimodel.github.io/EpiModel/reference/control.net.md).
  Key settings include `type` (disease model: `"SI"`, `"SIR"`, or
  `"SIS"`), `nsteps` (number of time steps), `nsims` (number of
  simulations), `tergmLite` (lightweight mode for performance), and
  `resimulate.network` (required for models with vital dynamics). For
  extended models, custom module functions are also passed here.

- scenarios_list:

  A list of scenarios to be run. Produced by the
  [`create_scenario_list()`](https://epimodel.github.io/EpiModel/reference/create_scenario_list.md)
  function

- n_rep:

  The number of replication to be run for each scenario.

- n_cores:

  The number of CPUs on which the simulations will be run.

- output_dir:

  The folder where the simulation files are to be stored.

- ...:

  for compatibility reasons

## Value

`NULL`, invisibly. The simulations are saved in `output_dir`.

## Details

The batches run one after the other in the current R session, while the
simulations of a batch run in parallel, one `netsim` call per simulation
(see
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)),
by default with `future::plan("multisession", workers = n_cores)`. Using
`control$future.use.plan <- future::tweak(<your plan>)` will bypass the
`multisession` setting.

The workers do not have the packages attached in the current session, so
module code written in a script must call package functions with
`pkg::`, such as
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).
Modules from a package find the functions it imports. See the Parallel
Simulation section of
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

The EpiModelHPC package provides `slurmworkflow` step templates that run
the same code on a SLURM cluster.

## Checkpointing

This function takes care of editing `.checkpoint.dir` to create unique
sub directories for each scenario, batch and simulation:
"\<.checkpoint.dir\>/sim\_\_\_\_\<batch_num\>/sim\_". The `control.net`
way of setting up checkpoints can be used transparently. Unless
`.checkpoint.keep` is `TRUE`, the batch directory is removed once the
batch file is saved.
