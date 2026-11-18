# Validate a Restart Pool

Checks that a restart pool (a directory made by
[`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md),
holding a `pool_index.csv`) can be used by
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md).
This function is not called by the simulation functions: call it after
making a pool, or before starting a long run from it.

## Usage

``` r
validate_restart_pool(path)
```

## Arguments

- path:

  Path to a restart pool directory, or to a single restart point file.

## Value

A `data.frame` with one row per pool element (invisibly): `pool_num`,
`file_path`, `n_runs`, `nsteps` (the number of time steps kept in the
element, `1` for a pool made by
[`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md))
and `epimodel_version`.

## Details

Errors are raised when:

- `path` is a directory without a `pool_index.csv`, or its
  `pool_index.csv` has no `pool_num` or `epimodel_version` column, or
  does not number the elements `1` to `N`;

- an element listed in `pool_index.csv` has no `<pool_num>.rds` file;

- an element is not a single simulation `netsim` restart point with the
  elements needed to restart (`param`, `control`, `nwparam`, `epi`,
  `run`, `coef.form` and `num.nw`);

- an element was made with another major.minor version of EpiModel than
  the one installed, according to `pool_index.csv` (checked by
  [`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)
  before running) or to the element itself (checked by `netsim`). Run
  this function with the EpiModel version the simulations will use (e.g.
  in the HPC `renv`);

- the elements do not all hold the same number of time steps: `netsim`
  requires `control$start` to be that number + 1 for each of them.

Warnings are raised when:

- the pool directory holds `<n>.rds` files not listed in
  `pool_index.csv`, which are never used;

- `path` is a single file holding several simulations:
  [`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)
  restarts every replicate from its first simulation.

## See also

[`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md),
[`netsim_path_wrapper()`](https://epimodel.github.io/EpiModel/reference/netsim_path_wrapper.md)
