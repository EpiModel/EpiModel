# Print Method for Network Model Simulations

Prints a detailed summary of a stochastic network model simulation
object, including simulation metadata, model parameters, output
variables, and (optionally) network formation, duration, and dissolution
statistics with target comparisons.

## Usage

``` r
# S3 method for class 'netsim'
print(x, nwstats = TRUE, digits = 3, network = 1, ...)
```

## Arguments

- x:

  An object of class `netsim`, from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- nwstats:

  If `TRUE` (the default), print network statistics tables (formation,
  duration, and dissolution) when available.

- digits:

  Number of digits to print in the network statistics tables.

- network:

  Integer index of the network for which to display statistics, for
  multi-network models. Default is `1`.

- ...:

  Additional arguments (currently ignored).

## Details

Given a `netsim` object, `print.netsim` displays the following sections:

**Simulation summary**: model class, model type (e.g., SI, SIR, SIS),
number of simulations, number of time steps, and number of network
groups.

**Model parameters**: printed via `print.param.net()`, showing fixed and
(if applicable) random parameters.

**Model functions**: for extension models (`type = NULL`), lists the
names of all custom module functions.

**Model output**: the names of all epidemic output variables, saved
network objects, transmission matrices, and any other saved elements.

**Network statistics** (when `nwstats = TRUE` and network statistics
were saved during the simulation):

The formation, duration, and dissolution statistics tables are computed
and displayed in the same way as for
[`print.netdx()`](http://epimodel.github.io/EpiModel/reference/print.netdx.md).
Each table contains the columns `Target`, `Sim Mean`, `Pct Diff`,
`Sim SE`, `Z Score`, `SD(Sim Means)`, and `SD(Statistic)`. The
`Sim Mean` column is the mean statistic value across all time steps and
simulations. `Sim SE` is the standard error estimated using
[coda::effectiveSize](https://rdrr.io/pkg/coda/man/effectiveSize.html).
`Pct Diff` gives `(Sim Mean - Target) / Target` and `Z Score` gives
`(Sim Mean - Target) / Sim SE` when a target is available.
`SD(Sim Means)` is the standard deviation of per-simulation means, and
`SD(Statistic)` is the overall standard deviation of the statistic.

*Formation statistics*: each row corresponds to a network statistic from
the formation formula (e.g., `edges`, `nodematch`), compared against the
target statistics from network estimation. These statistics assess
whether the network structure is maintained at target levels during the
epidemic simulation (which is particularly important in open-population
models where demographic turnover can shift network structure).

*Duration statistics*: each row corresponds to a dissolution model dyad
type. In a homogeneous dissolution model (`~offset(edges)`), all dyads
are of the `edges` type. The statistic of interest is the mean age of
extant edges, compared against the target duration from
[`dissolution_coefs()`](http://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).

*Dissolution statistics*: same row structure as the duration table. The
statistic of interest is the edge dissolution rate (proportion of edges
dissolving per time step), compared against `1 / target duration`.

Duration and dissolution tables are only available when
`control$save.diss.stats = TRUE`, `control$save.network = TRUE`,
`control$tergmLite = FALSE`, and the dissolution formula is
`~offset(edges)`. When these conditions are not met, a note listing the
requirements is printed instead.

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md),
[`print.netdx()`](http://epimodel.github.io/EpiModel/reference/print.netdx.md),
[`summary.netsim()`](http://epimodel.github.io/EpiModel/reference/summary.netsim.md).
