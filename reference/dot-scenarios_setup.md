# Set Up the Arguments to Run a Batch of Scenario Simulations

Exported for the EpiModelHPC `slurmworkflow` step templates, which run
the same code on a SLURM cluster; it is not meant to be called directly.

## Usage

``` r
.scenarios_setup(
  path_to_x,
  param,
  init,
  control,
  scenarios_list,
  n_rep,
  n_cores,
  output_dir
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

## Value

a list of arguments for `.run_batch`

## Details

Expands the scenarios over their batches and builds the `MoreArgs` list
passed to
[`.run_batch()`](https://epimodel.github.io/EpiModel/reference/dot-run_batch.md)
by
[`netsim_scenarios()`](https://epimodel.github.io/EpiModel/reference/netsim_scenarios.md),
checking the `control` settings the scenario runners support.
