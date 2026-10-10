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
#> 2    1    2   497     4 501       3       0
#> 3    1    3   495     6 501       3       1
#> 4    1    4   488    13 501      12       5
#> 5    1    5   468    33 501      21       1
#> 6    1    6   422    79 501      54       8
#> 7    1    7   354   147 501      92      24
#> 8    1    8   255   246 501     124      25
#> 9    1    9   157   344 501     148      50
#> 10   1   10   101   400 501     106      50
#> 11   2    1   500     1 501       0       0
#> 12   2    2   498     3 501       2       0
#> 13   2    3   496     5 501       2       0
#> 14   2    4   492     9 501       8       4
#> 15   2    5   483    18 501      12       3
#> 16   2    6   452    49 501      34       3
#> 17   2    7   382   119 501      82      12
#> 18   2    8   305   196 501     106      29
#> 19   2    9   221   280 501     128      44
#> 20   2   10   139   362 501     121      39
#> 21   3    1   500     1 501       0       0
#> 22   3    2   498     3 501       2       0
#> 23   3    3   495     6 501       4       1
#> 24   3    4   490    11 501       5       0
#> 25   3    5   478    23 501      13       1
#> 26   3    6   447    54 501      36       5
#> 27   3    7   386   115 501      70       9
#> 28   3    8   299   202 501     109      22
#> 29   3    9   189   312 501     139      29
#> 30   3   10   117   384 501     118      46
as.data.frame(mod, sim = 2)
#>    sim time s.num i.num num si.flow is.flow
#> 1    2    1   500     1 501       0       0
#> 2    2    2   498     3 501       2       0
#> 3    2    3   496     5 501       2       0
#> 4    2    4   492     9 501       8       4
#> 5    2    5   483    18 501      12       3
#> 6    2    6   452    49 501      34       3
#> 7    2    7   382   119 501      82      12
#> 8    2    8   305   196 501     106      29
#> 9    2    9   221   280 501     128      44
#> 10   2   10   139   362 501     121      39

# Time-specific means across simulations
as.data.frame(mod, out = "mean")
#>    time    s.num      i.num num    si.flow    is.flow
#> 1     1 500.0000   1.000000 501   0.000000  0.0000000
#> 2     2 497.6667   3.333333 501   2.333333  0.0000000
#> 3     3 495.3333   5.666667 501   3.000000  0.6666667
#> 4     4 490.0000  11.000000 501   8.333333  3.0000000
#> 5     5 476.3333  24.666667 501  15.333333  1.6666667
#> 6     6 440.3333  60.666667 501  41.333333  5.3333333
#> 7     7 374.0000 127.000000 501  81.333333 15.0000000
#> 8     8 286.3333 214.666667 501 113.000000 25.3333333
#> 9     9 189.0000 312.000000 501 138.333333 41.0000000
#> 10   10 119.0000 382.000000 501 115.000000 45.0000000

# Time-specific standard deviations across simulations
as.data.frame(mod, out = "sd")
#>    time      s.num      i.num num    si.flow    is.flow
#> 1     1  0.0000000  0.0000000   0  0.0000000  0.0000000
#> 2     2  0.5773503  0.5773503   0  0.5773503  0.0000000
#> 3     3  0.5773503  0.5773503   0  1.0000000  0.5773503
#> 4     4  2.0000000  2.0000000   0  3.5118846  2.6457513
#> 5     5  7.6376262  7.6376262   0  4.9328829  1.1547005
#> 6     6 16.0727513 16.0727513   0 11.0151411  2.5166115
#> 7     7 17.4355958 17.4355958   0 11.0151411  7.9372539
#> 8     8 27.3007936 27.3007936   0  9.6436508  3.5118846
#> 9     9 32.0000000 32.0000000   0 10.0166528 10.8166538
#> 10   10 19.0787840 19.0787840   0  7.9372539  5.5677644

# Time-specific quantile values across simulations
as.data.frame(mod, out = "qnt", qval = 0.25)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 497.5   3.0 501     2.0     0.0
#> 3     3 495.0   5.5 501     2.5     0.5
#> 4     4 489.0  10.0 501     6.5     2.0
#> 5     5 473.0  20.5 501    12.5     1.0
#> 6     6 434.5  51.5 501    35.0     4.0
#> 7     7 368.0 117.0 501    76.0    10.5
#> 8     8 277.0 199.0 501   107.5    23.5
#> 9     9 173.0 296.0 501   133.5    36.5
#> 10   10 109.0 373.0 501   112.0    42.5
as.data.frame(mod, out = "qnt", qval = 0.75)
#>    time s.num i.num num si.flow is.flow
#> 1     1 500.0   1.0 501     0.0     0.0
#> 2     2 498.0   3.5 501     2.5     0.0
#> 3     3 495.5   6.0 501     3.5     1.0
#> 4     4 491.0  12.0 501    10.0     4.5
#> 5     5 480.5  28.0 501    17.0     2.0
#> 6     6 449.5  66.5 501    45.0     6.5
#> 7     7 384.0 133.0 501    87.0    18.0
#> 8     8 302.0 224.0 501   116.5    27.0
#> 9     9 205.0 328.0 501   143.5    47.0
#> 10   10 128.0 392.0 501   119.5    48.0

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
