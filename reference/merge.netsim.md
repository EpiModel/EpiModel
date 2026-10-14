# Merge Model Simulations across netsim Objects

Merges epidemiological data from two independent simulations of
stochastic network models from `netsim`.

## Usage

``` r
# S3 method for class 'netsim'
merge(
  x,
  y,
  keep.transmat = TRUE,
  keep.network = TRUE,
  keep.nwstats = TRUE,
  keep.other = TRUE,
  param.error = TRUE,
  keep.diss.stats = TRUE,
  keep.run = TRUE,
  keep.cumulative.edgelist = FALSE,
  keep.attr.history = TRUE,
  keep.raw.records = TRUE,
  ...
)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- y:

  Another `EpiModel` object of class
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md),
  with the identical model parameterization as `x`.

- keep.transmat:

  If `TRUE`, keep the transmission matrices from the original `x` and
  `y` elements. Note: transmission matrices only saved when
  (`save.transmat == TRUE`).

- keep.network:

  If `TRUE`, keep the `networkDynamic` objects from the original `x` and
  `y` elements. Note: network only saved when (`tergmLite == FALSE`).

- keep.nwstats:

  If `TRUE`, keep the network statistics (as set by the
  `nwstats.formula` parameter in `control.net`) from the original `x`
  and `y` elements.

- keep.other:

  If `TRUE`, keep the other simulation elements (as set by the
  `save.other` parameter in `control.net`) from the original `x` and `y`
  elements. Elements with their own `keep.*` argument (e.g. `run`)
  follow that argument instead.

- param.error:

  If `TRUE`, if `x` and `y` have different params (in
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md))
  or controls (passed in
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md))
  an error will prevent the merge. Use `FALSE` to override that check.
  The controls for how the runs were executed rather than what they
  simulated (`nsims`, `ncores`, `verbose`, `verbose.int`, and the
  checkpoint settings) may differ.

- keep.diss.stats:

  If `TRUE`, keep `diss.stats` from the original `x` and `y` objects.

- keep.run:

  If `TRUE`, keep the `run` sublists (as set by the `save.run` parameter
  in `control.net`) from the original `x` and `y` elements. These are
  required to restart from a simulation of the merged object in
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md),
  selected with
  [`get_sims()`](https://epimodel.github.io/EpiModel/reference/get_sims.md).

- keep.cumulative.edgelist:

  If `TRUE`, keep the cumulative edgelists (as set by the
  `save.cumulative.edgelist` parameter in `control.net`) from the
  original `x` and `y` elements. `FALSE` by default, as these grow with
  the length of the simulation.

- keep.attr.history:

  If `TRUE`, keep the attribute histories (as set by
  [`record_attr_history()`](https://epimodel.github.io/EpiModel/reference/record_attr_history.md))
  from the original `x` and `y` elements.

- keep.raw.records:

  If `TRUE`, keep the raw records (as set by
  [`record_raw_object()`](https://epimodel.github.io/EpiModel/reference/record_raw_object.md))
  from the original `x` and `y` elements.

- ...:

  Additional merge arguments (not currently used).

## Value

An `EpiModel` object of class
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
containing the data from both `x` and `y`.

## Details

This merge function combines the results of two independent simulations
of [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
class models, simulated under separate function calls. The model
parameterization between the two calls must be exactly the same, except
for the number of simulations in each call and the execution settings
listed under `param.error`. This allows for manual parallelization of
model simulations.

This merge function does not work the same as the default merge, which
allows for a combined object where the structure differs between the
input elements. Instead, the function checks that objects are identical
in model parameterization in every respect (except number of
simulations) and binds the results.

After dropping the elements not kept (see the `keep.*` arguments), `x`
and `y` must hold the same elements.

## Examples

``` r
# \donttest{
# Network model
nw <- network_initialize(n = 100)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 10)
est <- netest(nw,
  formation = ~edges, target.stats = 25,
  coef.diss = coef.diss, verbose = FALSE
)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Epidemic models
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(
  type = "SI", nsteps = 20, nsims = 2,
  save.nwstats = TRUE,
  nwstats.formula = ~ edges + degree(0),
  verbose = FALSE
)
x <- netsim(est, param, init, control)
y <- netsim(est, param, init, control)

# Merging
z <- merge(x, y)

# Examine separate and merged data
as.data.frame(x)
#>    sim time s.num i.num num si.flow
#> 1    1    1    99     1 100      NA
#> 2    1    2    97     3 100       2
#> 3    1    3    95     5 100       2
#> 4    1    4    95     5 100       0
#> 5    1    5    95     5 100       0
#> 6    1    6    94     6 100       1
#> 7    1    7    94     6 100       0
#> 8    1    8    92     8 100       2
#> 9    1    9    92     8 100       0
#> 10   1   10    92     8 100       0
#> 11   1   11    91     9 100       1
#> 12   1   12    91     9 100       0
#> 13   1   13    91     9 100       0
#> 14   1   14    90    10 100       1
#> 15   1   15    89    11 100       1
#> 16   1   16    87    13 100       2
#> 17   1   17    86    14 100       1
#> 18   1   18    86    14 100       0
#> 19   1   19    86    14 100       0
#> 20   1   20    86    14 100       0
#> 21   2    1    99     1 100      NA
#> 22   2    2    98     2 100       1
#> 23   2    3    98     2 100       0
#> 24   2    4    98     2 100       0
#> 25   2    5    98     2 100       0
#> 26   2    6    98     2 100       0
#> 27   2    7    98     2 100       0
#> 28   2    8    98     2 100       0
#> 29   2    9    98     2 100       0
#> 30   2   10    98     2 100       0
#> 31   2   11    98     2 100       0
#> 32   2   12    98     2 100       0
#> 33   2   13    98     2 100       0
#> 34   2   14    98     2 100       0
#> 35   2   15    98     2 100       0
#> 36   2   16    97     3 100       1
#> 37   2   17    97     3 100       0
#> 38   2   18    97     3 100       0
#> 39   2   19    97     3 100       0
#> 40   2   20    97     3 100       0
as.data.frame(y)
#>    sim time s.num i.num num si.flow
#> 1    1    1    99     1 100      NA
#> 2    1    2    98     2 100       1
#> 3    1    3    98     2 100       0
#> 4    1    4    98     2 100       0
#> 5    1    5    98     2 100       0
#> 6    1    6    98     2 100       0
#> 7    1    7    98     2 100       0
#> 8    1    8    98     2 100       0
#> 9    1    9    98     2 100       0
#> 10   1   10    98     2 100       0
#> 11   1   11    98     2 100       0
#> 12   1   12    98     2 100       0
#> 13   1   13    98     2 100       0
#> 14   1   14    98     2 100       0
#> 15   1   15    98     2 100       0
#> 16   1   16    98     2 100       0
#> 17   1   17    98     2 100       0
#> 18   1   18    98     2 100       0
#> 19   1   19    98     2 100       0
#> 20   1   20    98     2 100       0
#> 21   2    1    99     1 100      NA
#> 22   2    2    99     1 100       0
#> 23   2    3    99     1 100       0
#> 24   2    4    99     1 100       0
#> 25   2    5    98     2 100       1
#> 26   2    6    98     2 100       0
#> 27   2    7    98     2 100       0
#> 28   2    8    98     2 100       0
#> 29   2    9    97     3 100       1
#> 30   2   10    97     3 100       0
#> 31   2   11    96     4 100       1
#> 32   2   12    96     4 100       0
#> 33   2   13    96     4 100       0
#> 34   2   14    96     4 100       0
#> 35   2   15    95     5 100       1
#> 36   2   16    95     5 100       0
#> 37   2   17    95     5 100       0
#> 38   2   18    95     5 100       0
#> 39   2   19    95     5 100       0
#> 40   2   20    95     5 100       0
as.data.frame(z)
#>    sim time s.num i.num num si.flow
#> 1    1    1    99     1 100      NA
#> 2    1    2    97     3 100       2
#> 3    1    3    95     5 100       2
#> 4    1    4    95     5 100       0
#> 5    1    5    95     5 100       0
#> 6    1    6    94     6 100       1
#> 7    1    7    94     6 100       0
#> 8    1    8    92     8 100       2
#> 9    1    9    92     8 100       0
#> 10   1   10    92     8 100       0
#> 11   1   11    91     9 100       1
#> 12   1   12    91     9 100       0
#> 13   1   13    91     9 100       0
#> 14   1   14    90    10 100       1
#> 15   1   15    89    11 100       1
#> 16   1   16    87    13 100       2
#> 17   1   17    86    14 100       1
#> 18   1   18    86    14 100       0
#> 19   1   19    86    14 100       0
#> 20   1   20    86    14 100       0
#> 21   2    1    99     1 100      NA
#> 22   2    2    98     2 100       1
#> 23   2    3    98     2 100       0
#> 24   2    4    98     2 100       0
#> 25   2    5    98     2 100       0
#> 26   2    6    98     2 100       0
#> 27   2    7    98     2 100       0
#> 28   2    8    98     2 100       0
#> 29   2    9    98     2 100       0
#> 30   2   10    98     2 100       0
#> 31   2   11    98     2 100       0
#> 32   2   12    98     2 100       0
#> 33   2   13    98     2 100       0
#> 34   2   14    98     2 100       0
#> 35   2   15    98     2 100       0
#> 36   2   16    97     3 100       1
#> 37   2   17    97     3 100       0
#> 38   2   18    97     3 100       0
#> 39   2   19    97     3 100       0
#> 40   2   20    97     3 100       0
#> 41   3    1    99     1 100      NA
#> 42   3    2    98     2 100       1
#> 43   3    3    98     2 100       0
#> 44   3    4    98     2 100       0
#> 45   3    5    98     2 100       0
#> 46   3    6    98     2 100       0
#> 47   3    7    98     2 100       0
#> 48   3    8    98     2 100       0
#> 49   3    9    98     2 100       0
#> 50   3   10    98     2 100       0
#> 51   3   11    98     2 100       0
#> 52   3   12    98     2 100       0
#> 53   3   13    98     2 100       0
#> 54   3   14    98     2 100       0
#> 55   3   15    98     2 100       0
#> 56   3   16    98     2 100       0
#> 57   3   17    98     2 100       0
#> 58   3   18    98     2 100       0
#> 59   3   19    98     2 100       0
#> 60   3   20    98     2 100       0
#> 61   4    1    99     1 100      NA
#> 62   4    2    99     1 100       0
#> 63   4    3    99     1 100       0
#> 64   4    4    99     1 100       0
#> 65   4    5    98     2 100       1
#> 66   4    6    98     2 100       0
#> 67   4    7    98     2 100       0
#> 68   4    8    98     2 100       0
#> 69   4    9    97     3 100       1
#> 70   4   10    97     3 100       0
#> 71   4   11    96     4 100       1
#> 72   4   12    96     4 100       0
#> 73   4   13    96     4 100       0
#> 74   4   14    96     4 100       0
#> 75   4   15    95     5 100       1
#> 76   4   16    95     5 100       0
#> 77   4   17    95     5 100       0
#> 78   4   18    95     5 100       0
#> 79   4   19    95     5 100       0
#> 80   4   20    95     5 100       0
# }
```
