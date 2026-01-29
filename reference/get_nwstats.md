# Extract Network Statistics from netsim or netdx Object

Extracts network statistics from a network epidemic model simulated with
`netsim` or a network diagnostics object simulated with `netdx`.
Statistics can be returned either as a single data frame or as a list of
matrices (one matrix for each simulation).

## Usage

``` r
get_nwstats(x, sim = NULL, network = 1, mode = c("data.frame", "list"))
```

## Arguments

- x:

  An `EpiModel` object of class
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
  or [`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md).

- sim:

  A vector of simulation numbers from the extracted object. (Default =
  NULL, all simulations are included)

- network:

  Network number, for `netsim` objects with multiple overlapping
  networks (advanced use, and not applicable to `netdx` objects).

- mode:

  Either `"data.frame"` or `"list"`, indicating the desired output.

## Value

A data frame or list of matrices containing the network statistics.

## Examples

``` r
# Two-group Bernoulli random graph TERGM
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

dx <- netdx(est, nsim = 3, nsteps = 10, verbose = FALSE,
            nwstats.formula = ~edges + isolates)
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range
get_nwstats(dx)
#>    time sim edges isolates
#> 1     1   1    38       50
#> 2     2   1    40       48
#> 3     3   1    41       48
#> 4     4   1    38       48
#> 5     5   1    40       47
#> 6     6   1    40       45
#> 7     7   1    39       47
#> 8     8   1    39       46
#> 9     9   1    39       45
#> 10   10   1    39       45
#> 11    1   2    44       41
#> 12    2   2    44       41
#> 13    3   2    46       36
#> 14    4   2    46       36
#> 15    5   2    46       36
#> 16    6   2    47       37
#> 17    7   2    50       36
#> 18    8   2    52       38
#> 19    9   2    57       37
#> 20   10   2    55       37
#> 21    1   3    57       31
#> 22    2   3    59       31
#> 23    3   3    60       28
#> 24    4   3    62       26
#> 25    5   3    63       25
#> 26    6   3    65       23
#> 27    7   3    65       21
#> 28    8   3    64       22
#> 29    9   3    62       24
#> 30   10   3    63       22
get_nwstats(dx, sim = 1)
#>    time sim edges isolates
#> 1     1   1    38       50
#> 2     2   1    40       48
#> 3     3   1    41       48
#> 4     4   1    38       48
#> 5     5   1    40       47
#> 6     6   1    40       45
#> 7     7   1    39       47
#> 8     8   1    39       46
#> 9     9   1    39       45
#> 10   10   1    39       45

# SI epidemic model
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3,
                       nwstats.formula = ~edges + meandeg + degree(0:5),
                       verbose = FALSE)
mod <- netsim(est, param, init, control)

# Extract the network statistics from all or sets of simulations
get_nwstats(mod)
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   1    56    1.12      35      31      22      11       1       0
#> 2     2   1    53    1.06      38      31      18      13       0       0
#> 3     3   1    57    1.14      36      33      14      15       2       0
#> 4     4   1    59    1.18      34      34      15      14       3       0
#> 5     5   1    59    1.18      34      36      11      16       3       0
#> 6     6   1    63    1.26      32      34      13      18       3       0
#> 7     7   1    63    1.26      33      30      18      16       3       0
#> 8     8   1    64    1.28      29      34      21      12       4       0
#> 9     9   1    64    1.28      26      37      23      11       3       0
#> 10   10   1    61    1.22      27      38      24       8       3       0
#> 11    1   2    61    1.22      35      30      21       9       3       1
#> 12    2   2    62    1.24      34      31      21       9       3       0
#> 13    3   2    56    1.12      38      32      17       9       2       1
#> 14    4   2    53    1.06      39      33      16      10       0       1
#> 15    5   2    52    1.04      39      31      20       8       1       1
#> 16    6   2    52    1.04      39      31      20       8       1       1
#> 17    7   2    54    1.08      38      30      22       7       2       1
#> 18    8   2    53    1.06      39      28      24       7       1       1
#> 19    9   2    54    1.08      38      28      26       6       0       2
#> 20   10   2    57    1.14      36      27      29       6       0       1
#> 21    1   3    50    1.00      31      44      20       4       1       0
#> 22    2   3    47    0.94      32      46      18       4       0       0
#> 23    3   3    47    0.94      31      48      17       4       0       0
#> 24    4   3    49    0.98      30      47      18       5       0       0
#> 25    5   3    47    0.94      34      44      16       6       0       0
#> 26    6   3    48    0.96      34      39      24       3       0       0
#> 27    7   3    47    0.94      39      33      23       5       0       0
#> 28    8   3    49    0.98      38      32      24       6       0       0
#> 29    9   3    51    1.02      36      32      27       4       1       0
#> 30   10   3    54    1.08      35      30      28       6       1       0
get_nwstats(mod, sim = 2)
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   2    61    1.22      35      30      21       9       3       1
#> 2     2   2    62    1.24      34      31      21       9       3       0
#> 3     3   2    56    1.12      38      32      17       9       2       1
#> 4     4   2    53    1.06      39      33      16      10       0       1
#> 5     5   2    52    1.04      39      31      20       8       1       1
#> 6     6   2    52    1.04      39      31      20       8       1       1
#> 7     7   2    54    1.08      38      30      22       7       2       1
#> 8     8   2    53    1.06      39      28      24       7       1       1
#> 9     9   2    54    1.08      38      28      26       6       0       2
#> 10   10   2    57    1.14      36      27      29       6       0       1
get_nwstats(mod, sim = c(1, 3))
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   1    56    1.12      35      31      22      11       1       0
#> 2     2   1    53    1.06      38      31      18      13       0       0
#> 3     3   1    57    1.14      36      33      14      15       2       0
#> 4     4   1    59    1.18      34      34      15      14       3       0
#> 5     5   1    59    1.18      34      36      11      16       3       0
#> 6     6   1    63    1.26      32      34      13      18       3       0
#> 7     7   1    63    1.26      33      30      18      16       3       0
#> 8     8   1    64    1.28      29      34      21      12       4       0
#> 9     9   1    64    1.28      26      37      23      11       3       0
#> 10   10   1    61    1.22      27      38      24       8       3       0
#> 11    1   3    50    1.00      31      44      20       4       1       0
#> 12    2   3    47    0.94      32      46      18       4       0       0
#> 13    3   3    47    0.94      31      48      17       4       0       0
#> 14    4   3    49    0.98      30      47      18       5       0       0
#> 15    5   3    47    0.94      34      44      16       6       0       0
#> 16    6   3    48    0.96      34      39      24       3       0       0
#> 17    7   3    47    0.94      39      33      23       5       0       0
#> 18    8   3    49    0.98      38      32      24       6       0       0
#> 19    9   3    51    1.02      36      32      27       4       1       0
#> 20   10   3    54    1.08      35      30      28       6       1       0

# On the fly summary stats
summary(get_nwstats(mod))
#>       time           sim        edges          meandeg         degree0     
#>  Min.   : 1.0   Min.   :1   Min.   :47.00   Min.   :0.940   Min.   :26.00  
#>  1st Qu.: 3.0   1st Qu.:1   1st Qu.:50.25   1st Qu.:1.005   1st Qu.:32.25  
#>  Median : 5.5   Median :2   Median :54.00   Median :1.080   Median :35.00  
#>  Mean   : 5.5   Mean   :2   Mean   :54.73   Mean   :1.095   Mean   :34.63  
#>  3rd Qu.: 8.0   3rd Qu.:3   3rd Qu.:59.00   3rd Qu.:1.180   3rd Qu.:38.00  
#>  Max.   :10.0   Max.   :3   Max.   :64.00   Max.   :1.280   Max.   :39.00  
#>     degree1         degree2         degree3          degree4     
#>  Min.   :27.00   Min.   :11.00   Min.   : 3.000   Min.   :0.000  
#>  1st Qu.:31.00   1st Qu.:17.25   1st Qu.: 6.000   1st Qu.:0.000  
#>  Median :32.50   Median :20.50   Median : 8.000   Median :1.000  
#>  Mean   :34.47   Mean   :20.33   Mean   : 8.667   Mean   :1.367  
#>  3rd Qu.:36.75   3rd Qu.:23.75   3rd Qu.:11.000   3rd Qu.:3.000  
#>  Max.   :48.00   Max.   :29.00   Max.   :18.000   Max.   :4.000  
#>     degree5      
#>  Min.   :0.0000  
#>  1st Qu.:0.0000  
#>  Median :0.0000  
#>  Mean   :0.3333  
#>  3rd Qu.:1.0000  
#>  Max.   :2.0000  
colMeans(get_nwstats(mod))
#>       time        sim      edges    meandeg    degree0    degree1    degree2 
#>  5.5000000  2.0000000 54.7333333  1.0946667 34.6333333 34.4666667 20.3333333 
#>    degree3    degree4    degree5 
#>  8.6666667  1.3666667  0.3333333 
```
