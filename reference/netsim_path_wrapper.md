# Run Single Simulations From a Restart Pool or a Single Object

Runs one
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
call with a single simulation for each element of `sim_nums`, in
parallel with
[`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
(using the current `future` plan).

## Usage

``` r
netsim_path_wrapper(path_to_x, param, init, control, sim_nums)
```

## Arguments

- path_to_x:

  Path to a restart pool directory (holding a `pool_index.csv`), or to a
  single file saved with `saveRDS` (a fitted network model or a restart
  point, see the `x` argument to
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)).

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

- sim_nums:

  The simulation numbers to run, used to pick the pool elements.

## Value

A list of `netsim` objects with one simulation each, in the order of
`sim_nums`.

## Details

If `path_to_x` is a restart pool directory (one holding a
`pool_index.csv`, see
[`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md))
of `N` elements, simulation `k` starts from pool element
`(k - 1) %% N + 1`: the pool is reused in order when there are more
simulations than pool elements; the path of the restart point used is
saved in `$restart_point`. Before any simulation starts, the elements
used are checked to have been made with the installed major.minor
version of EpiModel, according to `pool_index.csv`. `netsim` checks the
version of each restart point again when it starts. If `path_to_x` is a
file, every simulation starts from it (from its first simulation if it
holds several).

`control$nsims`, `control$ncores` and `control$future.use.plan` are
overridden. As `netsim` does for several simulations, `control$verbose`
is set to `FALSE` when `sim_nums` has more than one element. If
`control$.checkpoint.dir` is set, each simulation checkpoints in its own
`<.checkpoint.dir>/sim_<k>` sub directory.

Use
[`validate_restart_pool()`](https://epimodel.github.io/EpiModel/reference/validate_restart_pool.md)
to check a pool before using it.

## See also

[`make_restart_pool()`](https://epimodel.github.io/EpiModel/reference/make_restart_pool.md),
[`validate_restart_pool()`](https://epimodel.github.io/EpiModel/reference/validate_restart_pool.md)

## Examples

``` r
nw <- network_initialize(n = 50)
est <- netest(nw, formation = ~edges, target.stats = 20,
              coef.diss = dissolution_coefs(~offset(edges), 10, 0),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
est_path <- tempfile(fileext = ".rds")
saveRDS(est, est_path)

param <- param.net(inf.prob = 0.3, act.rate = 1)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 5, verbose = FALSE)

# Two simulations, each its own `netsim` call, on the current `future` plan
sims <- netsim_path_wrapper(est_path, param, init, control, sim_nums = 1:2)
length(sims)
#> [1] 2

unlink(est_path)
```
