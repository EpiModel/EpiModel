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
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  or
  [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md).

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
# \donttest{
# Two-group Bernoulli random graph TERGM
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

dx <- netdx(est, nsim = 3, nsteps = 10, verbose = FALSE,
            nwstats.formula = ~edges + isolates)
get_nwstats(dx)
#>    time sim edges isolates
#> 1     1   1    55       34
#> 2     2   1    53       34
#> 3     3   1    51       36
#> 4     4   1    50       39
#> 5     5   1    51       39
#> 6     6   1    51       40
#> 7     7   1    51       41
#> 8     8   1    51       40
#> 9     9   1    53       37
#> 10   10   1    54       36
#> 11    1   2    33       49
#> 12    2   2    36       45
#> 13    3   2    38       42
#> 14    4   2    36       45
#> 15    5   2    38       43
#> 16    6   2    37       48
#> 17    7   2    35       52
#> 18    8   2    40       50
#> 19    9   2    42       48
#> 20   10   2    39       52
#> 21    1   3    57       25
#> 22    2   3    54       29
#> 23    3   3    56       28
#> 24    4   3    58       26
#> 25    5   3    53       31
#> 26    6   3    55       31
#> 27    7   3    54       31
#> 28    8   3    52       31
#> 29    9   3    52       31
#> 30   10   3    56       28
get_nwstats(dx, sim = 1)
#>    time sim edges isolates
#> 1     1   1    55       34
#> 2     2   1    53       34
#> 3     3   1    51       36
#> 4     4   1    50       39
#> 5     5   1    51       39
#> 6     6   1    51       40
#> 7     7   1    51       41
#> 8     8   1    51       40
#> 9     9   1    53       37
#> 10   10   1    54       36

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
#> 1     1   1    55    1.10      33      37      20       7       3       0
#> 2     2   1    58    1.16      31      35      24       7       3       0
#> 3     3   1    58    1.16      31      34      26       6       3       0
#> 4     4   1    57    1.14      33      32      26       7       1       1
#> 5     5   1    57    1.14      33      32      26       7       1       1
#> 6     6   1    59    1.18      31      33      26       8       1       1
#> 7     7   1    61    1.22      28      36      24      10       2       0
#> 8     8   1    60    1.20      30      35      23      10       1       1
#> 9     9   1    60    1.20      29      36      23      11       0       1
#> 10   10   1    59    1.18      28      40      20      10       2       0
#> 11    1   2    48    0.96      39      36      18       5       1       1
#> 12    2   2    51    1.02      36      39      16       7       1       0
#> 13    3   2    47    0.94      38      40      14       6       2       0
#> 14    4   2    51    1.02      36      38      17       7       1       1
#> 15    5   2    50    1.00      34      41      17       7       1       0
#> 16    6   2    51    1.02      32      43      17       7       1       0
#> 17    7   2    47    0.94      34      46      13       6       1       0
#> 18    8   2    47    0.94      34      46      13       6       1       0
#> 19    9   2    54    1.08      28      47      15       9       1       0
#> 20   10   2    56    1.12      29      44      16       9       1       1
#> 21    1   3    43    0.86      43      35      16       5       1       0
#> 22    2   3    44    0.88      41      38      14       6       1       0
#> 23    3   3    41    0.82      44      37      12       7       0       0
#> 24    4   3    41    0.82      45      35      14       5       1       0
#> 25    5   3    42    0.84      45      33      16       5       1       0
#> 26    6   3    41    0.82      45      35      14       5       1       0
#> 27    7   3    42    0.84      46      34      12       6       2       0
#> 28    8   3    47    0.94      38      39      15       7       1       0
#> 29    9   3    43    0.86      40      42      11       6       1       0
#> 30   10   3    41    0.82      42      42       9       6       1       0
get_nwstats(mod, sim = 2)
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   2    48    0.96      39      36      18       5       1       1
#> 2     2   2    51    1.02      36      39      16       7       1       0
#> 3     3   2    47    0.94      38      40      14       6       2       0
#> 4     4   2    51    1.02      36      38      17       7       1       1
#> 5     5   2    50    1.00      34      41      17       7       1       0
#> 6     6   2    51    1.02      32      43      17       7       1       0
#> 7     7   2    47    0.94      34      46      13       6       1       0
#> 8     8   2    47    0.94      34      46      13       6       1       0
#> 9     9   2    54    1.08      28      47      15       9       1       0
#> 10   10   2    56    1.12      29      44      16       9       1       1
get_nwstats(mod, sim = c(1, 3))
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   1    55    1.10      33      37      20       7       3       0
#> 2     2   1    58    1.16      31      35      24       7       3       0
#> 3     3   1    58    1.16      31      34      26       6       3       0
#> 4     4   1    57    1.14      33      32      26       7       1       1
#> 5     5   1    57    1.14      33      32      26       7       1       1
#> 6     6   1    59    1.18      31      33      26       8       1       1
#> 7     7   1    61    1.22      28      36      24      10       2       0
#> 8     8   1    60    1.20      30      35      23      10       1       1
#> 9     9   1    60    1.20      29      36      23      11       0       1
#> 10   10   1    59    1.18      28      40      20      10       2       0
#> 11    1   3    43    0.86      43      35      16       5       1       0
#> 12    2   3    44    0.88      41      38      14       6       1       0
#> 13    3   3    41    0.82      44      37      12       7       0       0
#> 14    4   3    41    0.82      45      35      14       5       1       0
#> 15    5   3    42    0.84      45      33      16       5       1       0
#> 16    6   3    41    0.82      45      35      14       5       1       0
#> 17    7   3    42    0.84      46      34      12       6       2       0
#> 18    8   3    47    0.94      38      39      15       7       1       0
#> 19    9   3    43    0.86      40      42      11       6       1       0
#> 20   10   3    41    0.82      42      42       9       6       1       0

# On the fly summary stats
summary(get_nwstats(mod))
#>       time           sim        edges          meandeg         degree0     
#>  Min.   : 1.0   Min.   :1   Min.   :41.00   Min.   :0.820   Min.   :28.00  
#>  1st Qu.: 3.0   1st Qu.:1   1st Qu.:43.25   1st Qu.:0.865   1st Qu.:31.00  
#>  Median : 5.5   Median :2   Median :50.50   Median :1.010   Median :34.00  
#>  Mean   : 5.5   Mean   :2   Mean   :50.37   Mean   :1.007   Mean   :35.87  
#>  3rd Qu.: 8.0   3rd Qu.:3   3rd Qu.:57.00   3rd Qu.:1.140   3rd Qu.:40.75  
#>  Max.   :10.0   Max.   :3   Max.   :61.00   Max.   :1.220   Max.   :46.00  
#>     degree1         degree2         degree3      degree4         degree5      
#>  Min.   :32.00   Min.   : 9.00   Min.   : 5   Min.   :0.000   Min.   :0.0000  
#>  1st Qu.:35.00   1st Qu.:14.00   1st Qu.: 6   1st Qu.:1.000   1st Qu.:0.0000  
#>  Median :37.00   Median :16.00   Median : 7   Median :1.000   Median :0.0000  
#>  Mean   :38.00   Mean   :17.57   Mean   : 7   Mean   :1.267   Mean   :0.2667  
#>  3rd Qu.:40.75   3rd Qu.:22.25   3rd Qu.: 7   3rd Qu.:1.000   3rd Qu.:0.7500  
#>  Max.   :47.00   Max.   :26.00   Max.   :11   Max.   :3.000   Max.   :1.0000  
colMeans(get_nwstats(mod))
#>       time        sim      edges    meandeg    degree0    degree1    degree2 
#>  5.5000000  2.0000000 50.3666667  1.0073333 35.8666667 38.0000000 17.5666667 
#>    degree3    degree4    degree5 
#>  7.0000000  1.2666667  0.2666667 
# }
```
