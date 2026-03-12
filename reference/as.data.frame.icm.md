# Extract Model Data for Stochastic Models

This function extracts model simulations for objects of classes `icm`
and `netsim` into a data frame using the generic `as.data.frame`
function.

## Usage

``` r
# S3 method for class 'icm'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  out = "vals",
  sim = NULL,
  qval = NULL,
  repair = "drop",
  ...
)

# S3 method for class 'netsim'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  out = "vals",
  sim = NULL,
  repair = "drop",
  ...
)
```

## Arguments

- x:

  An `EpiModel` object of class `icm` or `netsim`.

- row.names:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

- optional:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

- out:

  Data output to data frame: `"mean"` for row means across simulations,
  `"sd"` for row standard deviations across simulations, `"qnt"` for row
  quantiles at the level specified in `qval`, or `"vals"` for values
  from individual simulations.

- sim:

  If `out="vals"`, the simulation number to output. If not specified,
  then data from all simulations will be output.

- qval:

  Quantile value required when `out="qnt"`.

- repair:

  What to do with epi trackers that are too short. "drop" will remove
  them from the output, "pad" will add `NA` rows to get to the right
  size(Default = "drop").

- ...:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

## Value

A data frame containing the data from `x`.

## Details

These methods work for both `icm` and `netsim` class models. The
available output includes time-specific means, standard deviations,
quantiles, and simulation values (compartment and flow sizes) from these
stochastic model classes. Means, standard deviations, and quantiles are
calculated by taking the row summary (i.e., each row of data corresponds
to a time step) across all simulations in the model output.

## Examples

``` r
## Stochastic ICM SIS model
param <- param.icm(inf.prob = 0.8, act.rate = 2, rec.rate = 0.1)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SIS", nsteps = 10,
                       nsims = 3, verbose = FALSE)
mod <- icm(param, init, control)

# Default output all simulation runs, default to all in stacked data.frame
as.data.frame(mod)
#>    sim time s.num i.num num si.flow is.flow
#> 1    1    1   500     1 501       0       0
#> 2    1    2   500     1 501       0       0
#> 3    1    3   500     1 501       1       1
#> 4    1    4   497     4 501       3       0
#> 5    1    5   493     8 501       7       3
#> 6    1    6   485    16 501      10       2
#> 7    1    7   462    39 501      29       6
#> 8    1    8   418    83 501      53       9
#> 9    1    9   326   175 501     112      20
#> 10   1   10   208   293 501     146      28
#> 11   2    1   500     1 501       0       0
#> 12   2    2   500     1 501       2       2
#> 13   2    3   500     1 501       0       0
#> 14   2    4   498     3 501       2       0
#> 15   2    5   493     8 501       6       1
#> 16   2    6   485    16 501       9       1
#> 17   2    7   470    31 501      21       6
#> 18   2    8   428    73 501      47       5
#> 19   2    9   351   150 501     101      24
#> 20   2   10   253   248 501     128      30
#> 21   3    1   500     1 501       0       0
#> 22   3    2   498     3 501       3       1
#> 23   3    3   495     6 501       4       1
#> 24   3    4   489    12 501       7       1
#> 25   3    5   471    30 501      22       4
#> 26   3    6   433    68 501      46       8
#> 27   3    7   354   147 501      98      19
#> 28   3    8   255   246 501     132      33
#> 29   3    9   150   351 501     141      36
#> 30   3   10    82   419 501     108      40
as.data.frame(mod, sim = 2)
#>    sim time s.num i.num num si.flow is.flow
#> 1    2    1   500     1 501       0       0
#> 2    2    2   500     1 501       2       2
#> 3    2    3   500     1 501       0       0
#> 4    2    4   498     3 501       2       0
#> 5    2    5   493     8 501       6       1
#> 6    2    6   485    16 501       9       1
#> 7    2    7   470    31 501      21       6
#> 8    2    8   428    73 501      47       5
#> 9    2    9   351   150 501     101      24
#> 10   2   10   253   248 501     128      30

# Time-specific means across simulations
as.data.frame(mod, out = "mean")
#>    time    s.num      i.num num    si.flow    is.flow
#> 1     1 500.0000   1.000000 501   0.000000  0.0000000
#> 2     2 499.3333   1.666667 501   1.666667  1.0000000
#> 3     3 498.3333   2.666667 501   1.666667  0.6666667
#> 4     4 494.6667   6.333333 501   4.000000  0.3333333
#> 5     5 485.6667  15.333333 501  11.666667  2.6666667
#> 6     6 467.6667  33.333333 501  21.666667  3.6666667
#> 7     7 428.6667  72.333333 501  49.333333 10.3333333
#> 8     8 367.0000 134.000000 501  77.333333 15.6666667
#> 9     9 275.6667 225.333333 501 118.000000 26.6666667
#> 10   10 181.0000 320.000000 501 127.333333 32.6666667

# Time-specific standard deviations across simulations
as.data.frame(mod, out = "sd")
#>    time      s.num      i.num num   si.flow    is.flow
#> 1     1   0.000000   0.000000   0  0.000000  0.0000000
#> 2     2   1.154701   1.154701   0  1.527525  1.0000000
#> 3     3   2.886751   2.886751   0  2.081666  0.5773503
#> 4     4   4.932883   4.932883   0  2.645751  0.5773503
#> 5     5  12.701706  12.701706   0  8.962886  1.5275252
#> 6     6  30.022214  30.022214   0 21.079216  3.7859389
#> 7     7  64.786830  64.786830   0 42.335958  7.5055535
#> 8     8  97.123633  97.123633   0 47.437678 15.1437556
#> 9     9 109.546033 109.546033   0 20.663978  8.3266640
#> 10   10  88.639720  88.639720   0 19.008770  6.4291005

# Time-specific quantile values across simulations
as.data.frame(mod, out = "qnt", qval = 0.25)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 499.0   1.0 501     1.0     0.5
#> 3     3 497.5   1.0 501     0.5     0.5
#> 4     4 493.0   3.5 501     2.5     0.0
#> 5     5 482.0   8.0 501     6.5     2.0
#> 6     6 459.0  16.0 501     9.5     1.5
#> 7     7 408.0  35.0 501    25.0     6.0
#> 8     8 336.5  78.0 501    50.0     7.0
#> 9     9 238.0 162.5 501   106.5    22.0
#> 10   10 145.0 270.5 501   118.0    29.0
as.data.frame(mod, out = "qnt", qval = 0.75)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 500.0   2.0 501     2.5     1.5
#> 3     3 500.0   3.5 501     2.5     1.0
#> 4     4 497.5   8.0 501     5.0     0.5
#> 5     5 493.0  19.0 501    14.5     3.5
#> 6     6 485.0  42.0 501    28.0     5.0
#> 7     7 466.0  93.0 501    63.5    12.5
#> 8     8 423.0 164.5 501    92.5    21.0
#> 9     9 338.5 263.0 501   126.5    30.0
#> 10   10 230.5 356.0 501   137.0    35.0

if (FALSE) { # \dontrun{
## Stochastic SI Network Model
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

param <- param.net(inf.prob = 0.5)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
mod <- netsim(est, param, init, control)

# Same data extraction methods as with ICMs
as.data.frame(mod)
as.data.frame(mod, sim = 2)
as.data.frame(mod, out = "mean")
as.data.frame(mod, out = "sd")
as.data.frame(mod, out = "qnt", qval = 0.25)
as.data.frame(mod, out = "qnt", qval = 0.75)
} # }
```
