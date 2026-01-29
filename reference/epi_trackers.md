# Function to run the user-provided epi trackers

see the "Working with Custom Attributes and Summary Statistics in
EpiModel" vignette.

## Usage

``` r
epi_trackers(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## The `tracker.list` list

`.tracker.list` is a list of NAMED functions stored in the `control`
list of the main `netsim_dat` class object.

## Tracker Functions

This function will apply the tracker functions present in the control
list `.tracker.list`. Each tracker must be a function with EXACTLY one
argument: the `netsim_dat` main list object. They must return a VALUE of
length one (numeric, logical or character).

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)

## Examples

``` r
if (FALSE) { # \dontrun{

# Create some trackers
epi_prop_infected <- function(dat) {
  # we need two attributes for our calculation: `status` and `active`
  needed_attributes <- c("status", "active")
  # we use `with` to simplify code
  output <- with(EpiModel::get_attr_list(dat, needed_attributes), {
    pop <- active == 1             # we only look at active nodes
    cond <- status == "i"   # which are infected
    # how many are `infected` among the `active`
    sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
  })
  return(output)
}

epi_s_num <- function(dat) {
  needed_attributes <- c("status")
  output <- with(get_attr_list(dat, needed_attributes), {
    sum(status == "s", na.rm = TRUE)
  })
  return(output)
}

# Store the trackers in a named list. The names will be used as column names
# for in the `epi` list
some.trackers <- list(
  prop_infected = epi_prop_infected,
  s_num         = epi_s_num
)

# Make a simple SI model with custom trackers
control <- EpiModel::control.net(
  type = "SI",
  nsims = 1,
  nsteps = 50,
  verbose = FALSE,
  .tracker.list = some.trackers
)

param <- EpiModel::param.net(
  inf.prob = 0.3,
  act.rate = 0.1
)

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- EpiModel::netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
  verbose = FALSE
)

init <- EpiModel::init.net(i.num = 10)
sim <- EpiModel::netsim(est, param, init, control)

d <- as.data.frame(sim)
d
} # }
```
