# Extract Degree Distributions from Network Diagnostics

Calculates the cumulative or momentary degree distribution from the
timed edgelist of a dynamic
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
simulation.

## Usage

``` r
get_degree_dist(
  x,
  type = c("cumulative", "momentary"),
  sims = NULL,
  window = NULL,
  unique.partners = TRUE
)
```

## Arguments

- x:

  An `EpiModel` object of class `netdx`, run with `dynamic = TRUE` and
  `keep.tedgelist = TRUE`.

- type:

  Distribution to calculate, either `"cumulative"` for the number of
  partners each node accumulates over the observation window, or
  `"momentary"` for the number of partners each node has at a single
  point in time, averaged over the time steps in the window.

- sims:

  A vector of simulation numbers to include. The default is all
  simulations in `x`.

- window:

  A vector of length 2 giving the first and last time step of the
  observation window, on the time scale of the timed edgelist (which
  starts at 0). The default is the full simulation, `c(0, x$nsteps)`.

- unique.partners:

  If `TRUE`, a dyad that forms, dissolves, and then re-forms within the
  window contributes one partner to the cumulative degree of each node;
  if `FALSE`, it contributes one per partnership. Only used when
  `type = "cumulative"`.

## Value

A `data.frame` with one row per simulation and degree value, and the
columns:

- `sim`: the simulation number.

- `degree`: the degree value, from 0 to the maximum observed across the
  selected simulations.

- `count`: the number of nodes with that degree. For
  `type = "momentary"`, this is the mean number of nodes per time step,
  and so need not be a whole number.

- `prop`: `count` divided by the network size.

## Details

The momentary degree distribution is the distribution plotted by
[`plot.netdx()`](https://epimodel.github.io/EpiModel/reference/plot.netdx.md)
when `degree(k)` terms are included in the `nwstats.formula` of
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md): it
counts the partners a node has at one moment in time. The cumulative
degree distribution counts the partners a node accumulates across the
whole observation window, so a node with one partner at a time and a
node with several simultaneous partners may have the same cumulative
degree.

Partnerships are counted if they are active at any point in the window:
a partnership contributes to the cumulative degree of both of its nodes
when its spell overlaps `window`, and to the momentary degree of both
nodes at each time step of the window in which it is active.
Partnerships that are ongoing when the simulation ends are counted for
the steps they are observed.

Because the cumulative degree grows with the length of the observation
window, `window` should be set to a substantively meaningful period (for
example, one year of time steps) when comparing across models or against
data.

## See also

[`plot.netdx()`](https://epimodel.github.io/EpiModel/reference/plot.netdx.md)
with `type = "cumldeg"` to plot these distributions, and
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) to
generate the input object.

## Examples

``` r
# \donttest{
nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 50,
              coef.diss = dissolution_coefs(~offset(edges), duration = 10),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
dx <- netdx(est, nsims = 2, nsteps = 50, keep.tedgelist = TRUE,
            verbose = FALSE)

# Cumulative degree over the full simulation
cd <- get_degree_dist(dx, type = "cumulative")
head(cd)
#>   sim degree count prop
#> 1   1      0     0 0.00
#> 2   1      1     2 0.02
#> 3   1      2     3 0.03
#> 4   1      3    10 0.10
#> 5   1      4    10 0.10
#> 6   1      5    18 0.18

# Mean cumulative degree per simulation
tapply(cd$degree * cd$prop, cd$sim, sum)
#>    1    2 
#> 5.88 5.66 

# Momentary degree over the last 25 steps
get_degree_dist(dx, type = "momentary", window = c(25, 50))
#>    sim degree       count         prop
#> 1    1      0 39.92307692 0.3992307692
#> 2    1      1 35.50000000 0.3550000000
#> 3    1      2 16.07692308 0.1607692308
#> 4    1      3  6.50000000 0.0650000000
#> 5    1      4  1.76923077 0.0176923077
#> 6    1      5  0.23076923 0.0023076923
#> 7    1      6  0.00000000 0.0000000000
#> 8    2      0 38.69230769 0.3869230769
#> 9    2      1 36.11538462 0.3611538462
#> 10   2      2 17.88461538 0.1788461538
#> 11   2      3  5.19230769 0.0519230769
#> 12   2      4  1.61538462 0.0161538462
#> 13   2      5  0.46153846 0.0046153846
#> 14   2      6  0.03846154 0.0003846154
# }
```
