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
  theme from the output, "pad" will add `NA` rows to get to the right
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
calculated by taking the row summary (i.e., each row of data is
corresponds to a time step) across all simulations in the model output.

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
#> 3    1    3   499     2 501       1       0
#> 4    1    4   494     7 501       5       0
#> 5    1    5   485    16 501      10       1
#> 6    1    6   465    36 501      22       2
#> 7    1    7   419    82 501      54       8
#> 8    1    8   342   159 501      89      12
#> 9    1    9   249   252 501     125      32
#> 10   1   10   139   362 501     146      36
#> 11   2    1   500     1 501       0       0
#> 12   2    2   499     2 501       1       0
#> 13   2    3   498     3 501       1       0
#> 14   2    4   496     5 501       3       1
#> 15   2    5   486    15 501      12       2
#> 16   2    6   460    41 501      30       4
#> 17   2    7   406    95 501      56       2
#> 18   2    8   303   198 501     120      17
#> 19   2    9   191   310 501     144      32
#> 20   2   10   110   391 501     119      38
#> 21   3    1   500     1 501       0       0
#> 22   3    2   501     0 501       0       1
#> 23   3    3   501     0 501       0       0
#> 24   3    4   501     0 501       0       0
#> 25   3    5   501     0 501       0       0
#> 26   3    6   501     0 501       0       0
#> 27   3    7   501     0 501       0       0
#> 28   3    8   501     0 501       0       0
#> 29   3    9   501     0 501       0       0
#> 30   3   10   501     0 501       0       0
as.data.frame(mod, sim = 2)
#>    sim time s.num i.num num si.flow is.flow
#> 1    2    1   500     1 501       0       0
#> 2    2    2   499     2 501       1       0
#> 3    2    3   498     3 501       1       0
#> 4    2    4   496     5 501       3       1
#> 5    2    5   486    15 501      12       2
#> 6    2    6   460    41 501      30       4
#> 7    2    7   406    95 501      56       2
#> 8    2    8   303   198 501     120      17
#> 9    2    9   191   310 501     144      32
#> 10   2   10   110   391 501     119      38

# Time-specific means across simulations
as.data.frame(mod, out = "mean")
#>    time    s.num      i.num num    si.flow    is.flow
#> 1     1 500.0000   1.000000 501  0.0000000  0.0000000
#> 2     2 500.0000   1.000000 501  0.3333333  0.3333333
#> 3     3 499.3333   1.666667 501  0.6666667  0.0000000
#> 4     4 497.0000   4.000000 501  2.6666667  0.3333333
#> 5     5 490.6667  10.333333 501  7.3333333  1.0000000
#> 6     6 475.3333  25.666667 501 17.3333333  2.0000000
#> 7     7 442.0000  59.000000 501 36.6666667  3.3333333
#> 8     8 382.0000 119.000000 501 69.6666667  9.6666667
#> 9     9 313.6667 187.333333 501 89.6666667 21.3333333
#> 10   10 250.0000 251.000000 501 88.3333333 24.6666667

# Time-specific standard deviations across simulations
as.data.frame(mod, out = "sd")
#>    time      s.num      i.num num    si.flow    is.flow
#> 1     1   0.000000   0.000000   0  0.0000000  0.0000000
#> 2     2   1.000000   1.000000   0  0.5773503  0.5773503
#> 3     3   1.527525   1.527525   0  0.5773503  0.0000000
#> 4     4   3.605551   3.605551   0  2.5166115  0.5773503
#> 5     5   8.962886   8.962886   0  6.4291005  1.0000000
#> 6     6  22.368132  22.368132   0 15.5349069  2.0000000
#> 7     7  51.507281  51.507281   0 31.7700068  4.1633320
#> 8     8 104.885652 104.885652   0 62.2923216  8.7368949
#> 9     9 164.806958 164.806958   0 78.2325593 18.4752086
#> 10   10 217.855457 217.855457   0 77.6809715 21.3853532

# Time-specific quantile values across simulations
as.data.frame(mod, out = "qnt", qval = 0.25)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 499.5   0.5 501     0.0     0.0
#> 3     3 498.5   1.0 501     0.5     0.0
#> 4     4 495.0   2.5 501     1.5     0.0
#> 5     5 485.5   7.5 501     5.0     0.5
#> 6     6 462.5  18.0 501    11.0     1.0
#> 7     7 412.5  41.0 501    27.0     1.0
#> 8     8 322.5  79.5 501    44.5     6.0
#> 9     9 220.0 126.0 501    62.5    16.0
#> 10   10 124.5 181.0 501    59.5    18.0
as.data.frame(mod, out = "qnt", qval = 0.75)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 500.5   1.5 501     0.5     0.5
#> 3     3 500.0   2.5 501     1.0     0.0
#> 4     4 498.5   6.0 501     4.0     0.5
#> 5     5 493.5  15.5 501    11.0     1.5
#> 6     6 483.0  38.5 501    26.0     3.0
#> 7     7 460.0  88.5 501    55.0     5.0
#> 8     8 421.5 178.5 501   104.5    14.5
#> 9     9 375.0 281.0 501   134.5    32.0
#> 10   10 320.0 376.5 501   132.5    37.0

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
