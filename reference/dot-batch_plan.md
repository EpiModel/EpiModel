# Build the `future` Plan to Run the Simulations of a Batch

Exported for the EpiModelHPC `slurmworkflow` step templates, which run
the same code on a SLURM cluster; it is not meant to be called directly.

## Usage

``` r
.batch_plan(control, strategy, n_cores)
```

## Arguments

- control:

  A `control.net` list. Its `future.use.plan` element is used as the
  plan when it holds one set with
  [`future::tweak()`](https://future.futureverse.org/reference/plan.html).

- strategy:

  A `future` backend function, e.g.
  [`future::multicore`](https://future.futureverse.org/reference/multicore.html).

- n_cores:

  The number of workers to use when `control$future.use.plan` holds no
  plan.

## Value

A `future` plan: `control$future.use.plan` when it is a
[`future::tweak()`](https://future.futureverse.org/reference/plan.html)
plan, else `future::tweak(strategy, workers = n_cores)`.

## Details

Only builds the plan: `with(future::plan(<plan>), local = TRUE)` must be
called by the function running the simulations, as it resets the plan
when its caller exits.

`strategy` is a function rather than a name, as
[`future::tweak()`](https://future.futureverse.org/reference/plan.html)
looks a name up from the calling frame (this namespace, then the search
path) before the `future` namespace.
