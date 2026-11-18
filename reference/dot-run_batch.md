# Run and Save One Batch of `netsim` Simulations

Exported for the EpiModelHPC `slurmworkflow` step templates, which run
the same code on a SLURM cluster; it is not meant to be called directly.

## Usage

``` r
.run_batch(
  scenario,
  batch_num,
  path_to_x,
  param,
  init,
  control,
  output_dir,
  n_batch,
  n_rep,
  n_cores
)
```

## Arguments

- scenario:

  A single "`EpiModel` scenario" to be used in the simulation

- batch_num:

  The batch number, calculated from the number of replications and CPUs
  required.

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

- output_dir:

  The folder where the simulation files are to be stored.

- n_batch:

  The number of batches to be run `ceiling(n_rep / n_cores)`.

- n_rep:

  The number of replication to be run for each scenario.

- n_cores:

  The number of CPUs on which the simulations will be run.

## Details

This inner function is called by
[`netsim_scenarios()`](https://epimodel.github.io/EpiModel/reference/netsim_scenarios.md).
It runs the simulations of batch `batch_num` with
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md),
one `netsim` call per simulation, merges them and saves them as
"sim\_\_\_\_\<batch_num\>.rds" in `output_dir`. The simulations of batch
`b` are numbered `(b - 1) * n_cores + 1` to `b * n_cores` over all the
batches of a scenario. With a restart pool of size `N`, simulation `k`
starts from pool element `(k - 1) %% N + 1`, the same for every
scenario. The simulations are run with the
[`future::tweak()`](https://future.futureverse.org/reference/plan.html)
plan in `control$future.use.plan`, or with
`future::plan("multisession", workers = n_cores)` if there is none. The
`slurmworkflow` step templates set it to `multicore` before calling this
function.

## Checkpointing

This function takes care of editing `.checkpoint.dir` to create unique
sub directories for each scenario, batch and simulation:
"\<.checkpoint.dir\>/sim\_\_\_\_\<batch_num\>/sim\_". The `control.net`
way of setting up checkpoints can be used transparently. Unless
`.checkpoint.keep` is `TRUE`, the batch directory is removed once the
batch file is saved.
