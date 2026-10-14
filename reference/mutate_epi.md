# Add New Epidemiology Variables

Inspired by
[`dplyr::mutate`](https://dplyr.tidyverse.org/reference/mutate.html),
`mutate_epi` adds new variables to the epidemiological and related
variables within simulated model objects of any class in `EpiModel`.

## Usage

``` r
mutate_epi(x, ...)
```

## Arguments

- x:

  An `EpiModel` object of class `dcm`, `icm`, or `netsim`.

- ...:

  Name-value pairs of expressions (see examples below).

## Value

The updated `EpiModel` object of class `dcm`, `icm`, or `netsim`.

## Examples

``` r
# DCM example
param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 500)
mod1 <- dcm(param, init, control)
mod1 <- mutate_epi(mod1, prev = i.num/num)
plot(mod1, y = "prev")


# \donttest{
# Network model example
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 1, i.num.g2 = 0)
control <- control.net(type = "SI", nsteps = 10, nsims = 3,
                       verbose = FALSE)
mod1 <- netsim(est1, param, init, control)
mod1
#> EpiModel Simulation
#> =======================
#> Model class: netsim
#> 
#> Simulation Summary
#> -----------------------
#> Model type: SI
#> No. simulations: 3
#> No. time steps: 10
#> No. NW groups: 2
#> 
#> Model Parameters
#> ---------------------------
#> inf.prob = 0.3
#> inf.prob.g2 = 0.15
#> act.rate = 1
#> groups = 2
#> 
#> Model Output
#> -----------------------
#> Variables: s.num i.num num s.num.g2 i.num.g2 num.g2 
#> si.flow si.flow.g2
#> Networks: sim1 ... sim3
#> Transmissions: sim1 ... sim3
#> 
#> Formation Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     50     49.5       -1  1.252  -0.399          6.52           5.9
#> 
#> 
#> Duration Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     20   17.318  -13.411  0.473  -5.671         1.364         1.798
#> 
#> Dissolution Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges   0.05    0.054    8.866  0.006    0.77         0.018         0.039
#> 

# Add the prevalences to the dataset
mod1 <- mutate_epi(mod1, i.prev = i.num / num,
                         i.prev.g2 = i.num.g2 / num.g2)
plot(mod1, y = c("i.prev", "i.prev.g2"), qnts = 0.5, legend = TRUE)


# Add incidence rate per 100 person years (assume time step = 1 week)
mod1 <- mutate_epi(mod1, ir100 = 5200*(si.flow + si.flow.g2) /
                                      (s.num + s.num.g2))
as.data.frame(mod1)
#>    sim time s.num i.num num s.num.g2 i.num.g2 num.g2 si.flow si.flow.g2 i.prev
#> 1    1    1    49     1  50       50        0     50      NA         NA   0.02
#> 2    1    2    48     2  50       50        0     50       1          0   0.04
#> 3    1    3    46     4  50       50        0     50       2          0   0.08
#> 4    1    4    45     5  50       50        0     50       1          0   0.10
#> 5    1    5    44     6  50       50        0     50       1          0   0.12
#> 6    1    6    44     6  50       49        1     50       0          1   0.12
#> 7    1    7    44     6  50       48        2     50       0          1   0.12
#> 8    1    8    43     7  50       48        2     50       1          0   0.14
#> 9    1    9    42     8  50       48        2     50       1          0   0.16
#> 10   1   10    41     9  50       48        2     50       1          0   0.18
#> 11   2    1    49     1  50       50        0     50      NA         NA   0.02
#> 12   2    2    49     1  50       50        0     50       0          0   0.02
#> 13   2    3    49     1  50       50        0     50       0          0   0.02
#> 14   2    4    49     1  50       50        0     50       0          0   0.02
#> 15   2    5    49     1  50       50        0     50       0          0   0.02
#> 16   2    6    49     1  50       50        0     50       0          0   0.02
#> 17   2    7    49     1  50       50        0     50       0          0   0.02
#> 18   2    8    49     1  50       50        0     50       0          0   0.02
#> 19   2    9    49     1  50       50        0     50       0          0   0.02
#> 20   2   10    49     1  50       50        0     50       0          0   0.02
#> 21   3    1    49     1  50       50        0     50      NA         NA   0.02
#> 22   3    2    49     1  50       50        0     50       0          0   0.02
#> 23   3    3    49     1  50       50        0     50       0          0   0.02
#> 24   3    4    48     2  50       50        0     50       1          0   0.04
#> 25   3    5    48     2  50       50        0     50       0          0   0.04
#> 26   3    6    48     2  50       50        0     50       0          0   0.04
#> 27   3    7    48     2  50       50        0     50       0          0   0.04
#> 28   3    8    48     2  50       50        0     50       0          0   0.04
#> 29   3    9    48     2  50       50        0     50       0          0   0.04
#> 30   3   10    48     2  50       50        0     50       0          0   0.04
#>    i.prev.g2     ir100
#> 1       0.00        NA
#> 2       0.00  53.06122
#> 3       0.00 108.33333
#> 4       0.00  54.73684
#> 5       0.00  55.31915
#> 6       0.02  55.91398
#> 7       0.04  56.52174
#> 8       0.04  57.14286
#> 9       0.04  57.77778
#> 10      0.04  58.42697
#> 11      0.00        NA
#> 12      0.00   0.00000
#> 13      0.00   0.00000
#> 14      0.00   0.00000
#> 15      0.00   0.00000
#> 16      0.00   0.00000
#> 17      0.00   0.00000
#> 18      0.00   0.00000
#> 19      0.00   0.00000
#> 20      0.00   0.00000
#> 21      0.00        NA
#> 22      0.00   0.00000
#> 23      0.00   0.00000
#> 24      0.00  53.06122
#> 25      0.00   0.00000
#> 26      0.00   0.00000
#> 27      0.00   0.00000
#> 28      0.00   0.00000
#> 29      0.00   0.00000
#> 30      0.00   0.00000
as.data.frame(mod1, out = "mean")
#>    time    s.num    i.num num s.num.g2  i.num.g2 num.g2   si.flow si.flow.g2
#> 1     1 49.00000 1.000000  50 50.00000 0.0000000     50       NaN        NaN
#> 2     2 48.66667 1.333333  50 50.00000 0.0000000     50 0.3333333  0.0000000
#> 3     3 48.00000 2.000000  50 50.00000 0.0000000     50 0.6666667  0.0000000
#> 4     4 47.33333 2.666667  50 50.00000 0.0000000     50 0.6666667  0.0000000
#> 5     5 47.00000 3.000000  50 50.00000 0.0000000     50 0.3333333  0.0000000
#> 6     6 47.00000 3.000000  50 49.66667 0.3333333     50 0.0000000  0.3333333
#> 7     7 47.00000 3.000000  50 49.33333 0.6666667     50 0.0000000  0.3333333
#> 8     8 46.66667 3.333333  50 49.33333 0.6666667     50 0.3333333  0.0000000
#> 9     9 46.33333 3.666667  50 49.33333 0.6666667     50 0.3333333  0.0000000
#> 10   10 46.00000 4.000000  50 49.33333 0.6666667     50 0.3333333  0.0000000
#>        i.prev   i.prev.g2    ir100
#> 1  0.02000000 0.000000000      NaN
#> 2  0.02666667 0.000000000 17.68707
#> 3  0.04000000 0.000000000 36.11111
#> 4  0.05333333 0.000000000 35.93269
#> 5  0.06000000 0.000000000 18.43972
#> 6  0.06000000 0.006666667 18.63799
#> 7  0.06000000 0.013333333 18.84058
#> 8  0.06666667 0.013333333 19.04762
#> 9  0.07333333 0.013333333 19.25926
#> 10 0.08000000 0.013333333 19.47566
# }
```
