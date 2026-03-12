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
#> 1     1   1    44       45
#> 2     2   1    46       42
#> 3     3   1    49       40
#> 4     4   1    51       38
#> 5     5   1    48       39
#> 6     6   1    46       42
#> 7     7   1    48       39
#> 8     8   1    49       39
#> 9     9   1    51       39
#> 10   10   1    52       39
#> 11    1   2    50       38
#> 12    2   2    51       35
#> 13    3   2    49       37
#> 14    4   2    50       37
#> 15    5   2    48       37
#> 16    6   2    49       38
#> 17    7   2    48       38
#> 18    8   2    49       38
#> 19    9   2    48       38
#> 20   10   2    49       37
#> 21    1   3    54       26
#> 22    2   3    55       27
#> 23    3   3    57       25
#> 24    4   3    53       27
#> 25    5   3    55       28
#> 26    6   3    56       27
#> 27    7   3    58       27
#> 28    8   3    57       29
#> 29    9   3    58       28
#> 30   10   3    59       28
get_nwstats(dx, sim = 1)
#>    time sim edges isolates
#> 1     1   1    44       45
#> 2     2   1    46       42
#> 3     3   1    49       40
#> 4     4   1    51       38
#> 5     5   1    48       39
#> 6     6   1    46       42
#> 7     7   1    48       39
#> 8     8   1    49       39
#> 9     9   1    51       39
#> 10   10   1    52       39

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
#> 1     1   1    47    0.94      40      33      22       4       0       1
#> 2     2   1    48    0.96      38      33      26       2       0       1
#> 3     3   1    46    0.92      38      37      22       2       0       1
#> 4     4   1    46    0.92      34      43      20       3       0       0
#> 5     5   1    48    0.96      34      41      20       5       0       0
#> 6     6   1    49    0.98      35      37      24       3       1       0
#> 7     7   1    51    1.02      33      38      24       4       1       0
#> 8     8   1    51    1.02      33      39      22       5       1       0
#> 9     9   1    50    1.00      32      42      21       4       1       0
#> 10   10   1    48    0.96      32      45      19       3       1       0
#> 11    1   2    52    1.04      30      43      22       3       2       0
#> 12    2   2    53    1.06      30      43      20       5       2       0
#> 13    3   2    51    1.02      32      43      17       7       1       0
#> 14    4   2    52    1.04      32      42      18       6       2       0
#> 15    5   2    51    1.02      33      41      18       7       1       0
#> 16    6   2    50    1.00      35      39      18       7       1       0
#> 17    7   2    50    1.00      36      37      19       7       1       0
#> 18    8   2    51    1.02      36      35      21       7       1       0
#> 19    9   2    54    1.08      32      37      23       7       1       0
#> 20   10   2    55    1.10      31      38      22       8       1       0
#> 21    1   3    50    1.00      35      40      17       6       2       0
#> 22    2   3    50    1.00      35      39      18       7       1       0
#> 23    3   3    51    1.02      34      40      17       8       1       0
#> 24    4   3    52    1.04      33      43      13       9       2       0
#> 25    5   3    54    1.08      33      39      17       9       2       0
#> 26    6   3    54    1.08      35      37      17       7       4       0
#> 27    7   3    54    1.08      34      41      14       5       6       0
#> 28    8   3    51    1.02      37      39      13       7       4       0
#> 29    9   3    46    0.92      42      35      15       5       3       0
#> 30   10   3    50    1.00      38      38      14       6       4       0
get_nwstats(mod, sim = 2)
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   2    52    1.04      30      43      22       3       2       0
#> 2     2   2    53    1.06      30      43      20       5       2       0
#> 3     3   2    51    1.02      32      43      17       7       1       0
#> 4     4   2    52    1.04      32      42      18       6       2       0
#> 5     5   2    51    1.02      33      41      18       7       1       0
#> 6     6   2    50    1.00      35      39      18       7       1       0
#> 7     7   2    50    1.00      36      37      19       7       1       0
#> 8     8   2    51    1.02      36      35      21       7       1       0
#> 9     9   2    54    1.08      32      37      23       7       1       0
#> 10   10   2    55    1.10      31      38      22       8       1       0
get_nwstats(mod, sim = c(1, 3))
#>    time sim edges meandeg degree0 degree1 degree2 degree3 degree4 degree5
#> 1     1   1    47    0.94      40      33      22       4       0       1
#> 2     2   1    48    0.96      38      33      26       2       0       1
#> 3     3   1    46    0.92      38      37      22       2       0       1
#> 4     4   1    46    0.92      34      43      20       3       0       0
#> 5     5   1    48    0.96      34      41      20       5       0       0
#> 6     6   1    49    0.98      35      37      24       3       1       0
#> 7     7   1    51    1.02      33      38      24       4       1       0
#> 8     8   1    51    1.02      33      39      22       5       1       0
#> 9     9   1    50    1.00      32      42      21       4       1       0
#> 10   10   1    48    0.96      32      45      19       3       1       0
#> 11    1   3    50    1.00      35      40      17       6       2       0
#> 12    2   3    50    1.00      35      39      18       7       1       0
#> 13    3   3    51    1.02      34      40      17       8       1       0
#> 14    4   3    52    1.04      33      43      13       9       2       0
#> 15    5   3    54    1.08      33      39      17       9       2       0
#> 16    6   3    54    1.08      35      37      17       7       4       0
#> 17    7   3    54    1.08      34      41      14       5       6       0
#> 18    8   3    51    1.02      37      39      13       7       4       0
#> 19    9   3    46    0.92      42      35      15       5       3       0
#> 20   10   3    50    1.00      38      38      14       6       4       0

# On the fly summary stats
summary(get_nwstats(mod))
#>       time           sim        edges          meandeg         degree0     
#>  Min.   : 1.0   Min.   :1   Min.   :46.00   Min.   :0.920   Min.   :30.00  
#>  1st Qu.: 3.0   1st Qu.:1   1st Qu.:49.25   1st Qu.:0.985   1st Qu.:32.25  
#>  Median : 5.5   Median :2   Median :51.00   Median :1.020   Median :34.00  
#>  Mean   : 5.5   Mean   :2   Mean   :50.50   Mean   :1.010   Mean   :34.40  
#>  3rd Qu.: 8.0   3rd Qu.:3   3rd Qu.:52.00   3rd Qu.:1.040   3rd Qu.:35.75  
#>  Max.   :10.0   Max.   :3   Max.   :55.00   Max.   :1.100   Max.   :42.00  
#>     degree1         degree2        degree3       degree4         degree5   
#>  Min.   :33.00   Min.   :13.0   Min.   :2.0   Min.   :0.000   Min.   :0.0  
#>  1st Qu.:37.00   1st Qu.:17.0   1st Qu.:4.0   1st Qu.:1.000   1st Qu.:0.0  
#>  Median :39.00   Median :19.0   Median :6.0   Median :1.000   Median :0.0  
#>  Mean   :39.23   Mean   :19.1   Mean   :5.6   Mean   :1.567   Mean   :0.1  
#>  3rd Qu.:41.75   3rd Qu.:22.0   3rd Qu.:7.0   3rd Qu.:2.000   3rd Qu.:0.0  
#>  Max.   :45.00   Max.   :26.0   Max.   :9.0   Max.   :6.000   Max.   :1.0  
colMeans(get_nwstats(mod))
#>      time       sim     edges   meandeg   degree0   degree1   degree2   degree3 
#>  5.500000  2.000000 50.500000  1.010000 34.400000 39.233333 19.100000  5.600000 
#>   degree4   degree5 
#>  1.566667  0.100000 
```
