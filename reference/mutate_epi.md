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


# Network model example
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
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
#> Fixed Parameters
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
#> edges     50   53.633    7.267  1.195   3.041         6.256         5.828
#> 
#> 
#> Duration Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     20   25.965   29.825  0.288  20.702         0.582         1.132
#> 
#> Dissolution Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges   0.05    0.037  -25.174  0.006   -2.14         0.006         0.032
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
#> 2    1    2    49     1  50       50        0     50       0          0   0.02
#> 3    1    3    49     1  50       50        0     50       0          0   0.02
#> 4    1    4    49     1  50       50        0     50       0          0   0.02
#> 5    1    5    49     1  50       50        0     50       0          0   0.02
#> 6    1    6    49     1  50       50        0     50       0          0   0.02
#> 7    1    7    49     1  50       50        0     50       0          0   0.02
#> 8    1    8    49     1  50       50        0     50       0          0   0.02
#> 9    1    9    49     1  50       50        0     50       0          0   0.02
#> 10   1   10    49     1  50       49        1     50       0          1   0.02
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
#> 24   3    4    49     1  50       50        0     50       0          0   0.02
#> 25   3    5    49     1  50       50        0     50       0          0   0.02
#> 26   3    6    49     1  50       49        1     50       0          1   0.02
#> 27   3    7    49     1  50       49        1     50       0          0   0.02
#> 28   3    8    49     1  50       49        1     50       0          0   0.02
#> 29   3    9    49     1  50       49        1     50       0          0   0.02
#> 30   3   10    49     1  50       49        1     50       0          0   0.02
#>    i.prev.g2    ir100
#> 1       0.00       NA
#> 2       0.00  0.00000
#> 3       0.00  0.00000
#> 4       0.00  0.00000
#> 5       0.00  0.00000
#> 6       0.00  0.00000
#> 7       0.00  0.00000
#> 8       0.00  0.00000
#> 9       0.00  0.00000
#> 10      0.02 53.06122
#> 11      0.00       NA
#> 12      0.00  0.00000
#> 13      0.00  0.00000
#> 14      0.00  0.00000
#> 15      0.00  0.00000
#> 16      0.00  0.00000
#> 17      0.00  0.00000
#> 18      0.00  0.00000
#> 19      0.00  0.00000
#> 20      0.00  0.00000
#> 21      0.00       NA
#> 22      0.00  0.00000
#> 23      0.00  0.00000
#> 24      0.00  0.00000
#> 25      0.00  0.00000
#> 26      0.02 53.06122
#> 27      0.02  0.00000
#> 28      0.02  0.00000
#> 29      0.02  0.00000
#> 30      0.02  0.00000
as.data.frame(mod1, out = "mean")
#>    time s.num i.num num s.num.g2  i.num.g2 num.g2 si.flow si.flow.g2 i.prev
#> 1     1    49     1  50 50.00000 0.0000000     50     NaN        NaN   0.02
#> 2     2    49     1  50 50.00000 0.0000000     50       0  0.0000000   0.02
#> 3     3    49     1  50 50.00000 0.0000000     50       0  0.0000000   0.02
#> 4     4    49     1  50 50.00000 0.0000000     50       0  0.0000000   0.02
#> 5     5    49     1  50 50.00000 0.0000000     50       0  0.0000000   0.02
#> 6     6    49     1  50 49.66667 0.3333333     50       0  0.3333333   0.02
#> 7     7    49     1  50 49.66667 0.3333333     50       0  0.0000000   0.02
#> 8     8    49     1  50 49.66667 0.3333333     50       0  0.0000000   0.02
#> 9     9    49     1  50 49.66667 0.3333333     50       0  0.0000000   0.02
#> 10   10    49     1  50 49.33333 0.6666667     50       0  0.3333333   0.02
#>      i.prev.g2    ir100
#> 1  0.000000000      NaN
#> 2  0.000000000  0.00000
#> 3  0.000000000  0.00000
#> 4  0.000000000  0.00000
#> 5  0.000000000  0.00000
#> 6  0.006666667 17.68707
#> 7  0.006666667  0.00000
#> 8  0.006666667  0.00000
#> 9  0.006666667  0.00000
#> 10 0.013333333 17.68707
```
