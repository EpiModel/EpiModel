# Dissolution Coefficients for Stochastic Network Models

Calculates dissolution coefficients, given a dissolution model and
average edge duration, to pass as offsets to an ERGM/TERGM model fit in
`netest`.

## Usage

``` r
dissolution_coefs(dissolution, duration, d.rate = 0)
```

## Arguments

- dissolution:

  Right-hand sided STERGM dissolution formula (see
  [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)).
  See below for list of supported dissolution models.

- duration:

  A vector of mean edge durations in arbitrary time units.

- d.rate:

  Departure or exit rate from the population, as a single homogeneous
  rate that applies to the entire population.

## Value

A list of class `disscoef` with the following elements:

- **dissolution:** right-hand sided STERGM dissolution formula passed in
  the function call.

- **duration:** mean edge durations passed into the function.

- **coef.crude:** mean durations transformed into logit coefficients.

- **coef.adj:** crude coefficients adjusted for the risk of departure on
  edge persistence, if the `d.rate` argument is supplied.

- **coef.form.corr:** corrections to be subtracted from formation
  coefficients.

- **d.rate:** the departure rate.

- **diss.model.type:** the form of the dissolution model; options
  include `edgesonly`, `nodematch`, and `nodemix`.

## Details

This function performs two calculations for dissolution coefficients
used in a network model estimated with
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md):

1.  **Transformation:** the mean durations of edges in a network are
    mathematically transformed to logit coefficients.

2.  **Adjustment:** in a dynamic network simulation in an open
    population (in which there are departures), it is further necessary
    to adjust these coefficients; this upward adjustment accounts for
    departure as a competing risk to edge dissolution.

The current dissolution models supported by this function and in network
model estimation in
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md) are
as follows:

- `~offset(edges)`: a homogeneous dissolution model in which the edge
  duration is the same for all partnerships. This requires specifying
  one duration value.

- `~offset(edges) + offset(nodematch("<attr>"))`: a heterogeneous model
  in which the edge duration varies by whether the nodes in the dyad
  have similar values of a specified attribute. The duration vector
  should now contain two values: the first is the mean edge duration of
  non-matched dyads, and the second is the duration of the matched
  dyads.

- `~offset(edges) + offset(nodemix("<attr>"))`: a heterogeneous model
  that extends the nodematch model to include non-binary attributes for
  homophily. The duration vector should first contain the base value,
  then the values for every other possible combination in the term.

## Examples

``` r
## Homogeneous dissolution model with no departures
dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#> Dissolution Coefficients
#> =======================
#> Dissolution Model: ~offset(edges)
#> Target Statistics: 25
#> Crude Coefficient: 3.178054
#> Mortality/Exit Rate: 0
#> Adjusted Coefficient: 3.178054

## Homogeneous dissolution model with departures
dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                  d.rate = 0.001)
#> Dissolution Coefficients
#> =======================
#> Dissolution Model: ~offset(edges)
#> Target Statistics: 25
#> Crude Coefficient: 3.178054
#> Mortality/Exit Rate: 0.001
#> Adjusted Coefficient: 3.229321

## Heterogeneous dissolution model in which same-race edges have
## shorter duration compared to mixed-race edges, with no departures
dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
                  duration = c(20, 10))
#> Dissolution Coefficients
#> =======================
#> Dissolution Model: ~offset(edges) + offset(nodematch("race"))
#> Target Statistics: 20 10
#> Crude Coefficient: 2.944439 -0.7472144
#> Mortality/Exit Rate: 0
#> Adjusted Coefficient: 2.944439 -0.7472144

## Heterogeneous dissolution model in which same-race edges have
## shorter duration compared to mixed-race edges, with departures
dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
                  duration = c(20, 10), d.rate = 0.001)
#> Dissolution Coefficients
#> =======================
#> Dissolution Model: ~offset(edges) + offset(nodematch("race"))
#> Target Statistics: 20 10
#> Crude Coefficient: 2.944439 -0.7472144
#> Mortality/Exit Rate: 0.001
#> Adjusted Coefficient: 2.98524 -0.7678231

if (FALSE) { # \dontrun{
## Extended example for differential homophily by age group
# Set up the network with nodes categorized into 5 age groups
nw <- network_initialize(n = 1000)
age.grp <- sample(1:5, 1000, TRUE)
nw <- set_vertex_attribute(nw, "age.grp", age.grp)

# durations = non-matched, age.grp1 & age.grp1, age.grp2 & age.grp2, ...
# TERGM will include differential homophily by age group with nodematch term
# Target stats for the formation model are overall edges, and then the number
# matched within age.grp 1, age.grp 2, ..., age.grp 5
form <- ~edges + nodematch("age.grp", diff = TRUE)
target.stats <- c(450, 100, 125, 40, 80, 100)

# Target stats for the dissolution model are duration of non-matched edges,
# then duration of edges matched within age.grp 1, age.grp 2, ..., age.grp 5
durs <- c(60, 30, 80, 100, 125, 160)
diss <- dissolution_coefs(~offset(edges) +
                            offset(nodematch("age.grp", diff = TRUE)),
                          duration = durs)

# Fit the TERGM
fit <- netest(nw, form, target.stats, diss)

# Full diagnostics to evaluate model fit
dx <- netdx(fit, nsims = 10, ncores = 4, nsteps = 300)
print(dx)

# Simulate one long time series to examine timed edgelist
dx <- netdx(fit, nsims = 1, nsteps = 5000, keep.tedgelist = TRUE)

# Extract timed-edgelist
te <- as.data.frame(dx)
head(te)

# Limit to non-censored edges
te <- te[which(te$onset.censored == FALSE & te$terminus.censored == FALSE),
         c("head", "tail", "duration")]
head(te)

# Look up the age group of head and tail nodes
te$ag.head <- age.grp[te$head]
te$ag.tail <- age.grp[te$tail]
head(te)

# Recover average edge durations for age-group pairing
mean(te$duration[te$ag.head != te$ag.tail])
mean(te$duration[te$ag.head == 1 & te$ag.tail == 1])
mean(te$duration[te$ag.head == 2 & te$ag.tail == 2])
mean(te$duration[te$ag.head == 3 & te$ag.tail == 3])
mean(te$duration[te$ag.head == 4 & te$ag.tail == 4])
mean(te$duration[te$ag.head == 5 & te$ag.tail == 5])
durs
} # }
```
