# Truncate Simulation Time Series

Left-truncates simulation epidemiological summary statistics and network
statistics at a specified time step.

## Usage

``` r
truncate_sim(x, at, reset.time)

# S3 method for class 'dcm'
truncate_sim(x, at, reset.time = TRUE)

# S3 method for class 'icm'
truncate_sim(x, at, reset.time = TRUE)

# S3 method for class 'netsim'
truncate_sim(x, at, reset.time = TRUE)
```

## Arguments

- x:

  Object of class `dcm`, `netsim`, or `icm`.

- at:

  Time step at which to left-truncate the time series.

- reset.time:

  If `TRUE`, the time step sequence in the truncated model will be reset
  to start at 1. If `FALSE`, the original time step values will be
  preserved. Default is `TRUE`.

## Value

The updated object of class `dcm`, `netsim`, or `icm`.

## Details

This function would be used when running a follow-up simulation from
time steps `b` to `c` after a burn-in period from time `a` to `b`, where
the final time window of interest for data analysis is `b` to `c` only.

For `icm` and `netsim` objects, `at` is a time step as labeled by
[`as.data.frame.icm()`](https://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md).
With `reset.time = FALSE`, the time step of the first row kept is stored
in `control$epi.start`, from which
[`as.data.frame.icm()`](https://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md),
the `plot` and `summary` methods, and
[`comp_plot()`](https://epimodel.github.io/EpiModel/reference/comp_plot.md)
label the time steps. Without it, the rows are labeled from 1, as in the
output of
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md) and
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md),
including restarted simulations.

For `netsim` objects, the network statistics are truncated along with
the epi trackers, as both hold one row per time step. With
`reset.time = TRUE`, the times of the attribute history, the raw
records, the transmission matrix, and the cumulative edgelist are
shifted to the new time steps as well, so that
[`get_attr_at()`](https://epimodel.github.io/EpiModel/reference/get_attr_at.md),
[`get_transmat()`](https://epimodel.github.io/EpiModel/reference/get_transmat.md)
and
[`make_networkDynamic()`](https://epimodel.github.io/EpiModel/reference/make_networkDynamic.md)
read them in the same frame. Their records before `at` are kept, with
times of 0 or less, as the state of the nodes and edges at `at` depends
on them and the transmissions before `at` remain part of the
transmission chains.

The saved networks and the `run` sublists, which hold the state needed
to restart a simulation, keep their original time steps. Build a restart
point with
[`make_restart_point()`](https://epimodel.github.io/EpiModel/reference/make_restart_point.md)
from the object before truncating it.

## Examples

``` r
# DCM examples
param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 200)
mod1 <- dcm(param, init, control)
plot(mod1)


# Reset time
mod2a <- truncate_sim(mod1, at = 150)
plot(mod2a)

head(as.data.frame(mod2a))
#>   time     s.num    i.num  si.flow num
#> 1    1 112.84477 388.1552 4.311259 501
#> 2    2 108.53351 392.4665 4.190824 501
#> 3    3 104.34268 396.6573 4.070347 501
#> 4    4 100.27233 400.7277 3.950123 501
#> 5    5  96.32221 404.6778 3.830429 501
#> 6    6  92.49178 408.5082 3.711522 501

# Do not reset time
mod2b <- truncate_sim(mod1, at = 150, reset.time = FALSE)
plot(mod2b)

head(as.data.frame(mod2b))
#>   time     s.num    i.num  si.flow num
#> 1  150 112.84477 388.1552 4.311259 501
#> 2  151 108.53351 392.4665 4.190824 501
#> 3  152 104.34268 396.6573 4.070347 501
#> 4  153 100.27233 400.7277 3.950123 501
#> 5  154  96.32221 404.6778 3.830429 501
#> 6  155  92.49178 408.5082 3.711522 501

# \donttest{
# ICM example
param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SI", nsteps = 200, nsims = 1)
mod1 <- icm(param, init, control)

# Reset time
mod2a <- truncate_sim(mod1, at = 150)
plot(mod2a)

head(as.data.frame(mod2a))
#>   sim time s.num i.num num si.flow
#> 1   1    1   109   392 501       4
#> 2   1    2   106   395 501       3
#> 3   1    3   101   400 501       5
#> 4   1    4    92   409 501       9
#> 5   1    5    89   412 501       3
#> 6   1    6    87   414 501       2

# Do not reset time
mod2b <- truncate_sim(mod1, at = 150, reset.time = FALSE)
plot(mod2b)

head(as.data.frame(mod2b))
#>   sim time s.num i.num num si.flow
#> 1   1  150   109   392 501       4
#> 2   1  151   106   395 501       3
#> 3   1  152   101   400 501       5
#> 4   1  153    92   409 501       9
#> 5   1  154    89   412 501       3
#> 6   1  155    87   414 501       2
# }
```
