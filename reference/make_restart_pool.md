# Make a Restart Pool From Simulations Run With Scenarios

Creates a directory of single simulation restart points (`1.rds`,
`2.rds`, ...) from the batch files written by
[`netsim_scenarios()`](https://epimodel.github.io/EpiModel/reference/netsim_scenarios.md).
Such a directory can be passed as `path_to_x` to the scenario runners
(see
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)).

## Usage

``` r
make_restart_pool(
  pool_dir,
  scenarios_dir,
  scenario_name,
  keep_sims,
  time_attrs
)
```

## Arguments

- pool_dir:

  Path to the restart pool directory to create. It must not exist.

- scenarios_dir:

  The directory containing the `sim__<scenario>__<batch>.rds` files.

- scenario_name:

  The name of the scenario to make the pool from.

- keep_sims:

  A `data.frame` with a `batch_number` and a `sim_number` column, one
  row per simulation to put in the pool. `sim_number` is the position of
  the simulation within its batch file (the `sim_number` column of
  [`merge_netsim_scenarios_tibble()`](https://epimodel.github.io/EpiModel/reference/merge_netsim_scenarios_tibble.md),
  not the `sim` one).

- time_attrs:

  A character vector of the attributes holding time steps, passed to
  [`make_restart_point()`](https://epimodel.github.io/EpiModel/reference/make_restart_point.md).

## Value

The pool index `data.frame` (invisibly).

## Details

Pool element `i` is made from row `i` of `keep_sims`, with
[`make_restart_point()`](https://epimodel.github.io/EpiModel/reference/make_restart_point.md)
(keeping only the last time step). Each batch file is read once, and
each element is written to the pool as soon as it is made, so only one
batch file's elements are in memory at a time. The pool directory also
receives a `pool_index.csv` file with one row per pool element:
`pool_num`, `scenario_name`, `batch_number`, `sim_number`, `file_path`
(the batch file it was made from) and `epimodel_version` (the EpiModel
version recorded in the element, which `netsim` requires to have the
installed major.minor version to restart from it). A directory is used
as a restart pool only if it holds this file.

The pool is written in a hidden `.<pool name>.partial` directory next to
`pool_dir`, renamed to `pool_dir` once complete: `pool_dir` never holds
a partial pool, even if R is killed. The partial directory is removed if
an error occurs, or at the next call if R was killed.

The EpiModelHPC package provides `slurmworkflow` step templates that run
[`netsim_scenarios()`](https://epimodel.github.io/EpiModel/reference/netsim_scenarios.md)
and this function on a SLURM cluster.

## See also

[`validate_restart_pool()`](https://epimodel.github.io/EpiModel/reference/validate_restart_pool.md),
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d_sim <- readRDS("data/run/calib/merged_tibbles/df__default.rds")
keep_sims <- dplyr::distinct(d_sim, batch_number, sim_number)[1:32, ]
make_restart_pool(
  pool_dir = "data/run/restart_pool",
  scenarios_dir = "data/run/calib",
  scenario_name = "default",
  keep_sims = keep_sims,
  time_attrs = c("infTime")
)
validate_restart_pool("data/run/restart_pool")
} # }
```
