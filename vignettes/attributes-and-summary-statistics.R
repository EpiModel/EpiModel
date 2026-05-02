## ----echo = FALSE, include = FALSE--------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE, include = FALSE-------------------------------------
library(EpiModel)

## -----------------------------------------------------------------------------
get_core_attributes()

## ----accessing-attr, eval = FALSE---------------------------------------------
# active <- get_attr(dat, "active")

## ----modyfing-attr, eval = FALSE----------------------------------------------
# aging_module <- function(dat, at) {
# 
#   # Extract current attribute
#   age <- get_attr(dat, "age")
# 
#   # Aging process
#   new_age <- age + 1
# 
#   # Output updated attributes
#   dat <- set_attr(dat, "age", new_age)
# 
#   return(dat)
# }

## ----recording_attr_hist, eval = FALSE----------------------------------------
# viral_load_logger_module <- function(dat, at) {
# 
#   # Run every 10 time steps
#   if (at %% 10 == 0) {
# 
#     # Attributes
#     status <- get_attr(dat, "status")
#     viral_load <- get_attr(dat, "viral_load")
# 
#     infected <- which(status == "i")
# 
#     dat <- record_attr_history(
#       dat, at,
#       "viral_load",
#       infected,
#       viral_load[infected]
#     )
#   }
# 
#   return(dat)
# }

## ----access_attr_hist, eval = FALSE-------------------------------------------
# sim <- netsim(est, param, init, control)
# attr_history <- get_attr_history(sim)

## ----get_attr_ex, eval = FALSE------------------------------------------------
# get_attr_history(sim)
# 
# # $viral_load
# #    sim step attribute    uids values
# # 1    1   10 viral_load   1001   2000
# # 2    1   10 viral_load   1002   1878
# # 3    1   20 viral_load   1001   1500
# # 4    1   20 viral_load   1002    300
# # ...
# #
# # $status
# #    sim step attribute  uids  values
# # 1    1   22 status     1001       i
# # 2    1   64 status     1002       i
# # 3    1  110 status     1001       r
# # 4    1  220 status     1002       r
# # ...

## ----epi-in-modules, eval = FALSE---------------------------------------------
# aging_track_module <- function(dat, at) {
# 
#   # Attributes
#   age <- get_attr(dat, "age")
# 
#   # Aging process
#   new_age <- age + 1
# 
#   # Calculate summary statistics
#   mean_age <- mean(new_age)
#   prev_mean_age <- get_epi(dat, at - 1, "mean_age")
#   age_change <- mean_age - prev_mean_age
# 
#   # Update nodal attributes
#   dat <- set_attr(dat, "age", new_age)
# 
#   # Update epidemic trackers
#   dat <- set_epi(dat, "mean_age", at, mean_age)
#   dat <- set_epi(dat, "age_change", at, age_change)
# 
#   return(dat)
# }

## ----eval = FALSE-------------------------------------------------------------
# sim <- netsim(est, param, init, control)
# 
# # Raw per-simulation data
# as.data.frame(sim)
# 
# # Means across simulations
# as.data.frame(sim, out = "mean")
# 
# # Plot and summary
# plot(sim)
# summary(sim, at = 50)

## ----eval = FALSE-------------------------------------------------------------
# # Add incidence rate and prevalence
# sim <- mutate_epi(sim, prev = i.num / num)
# sim <- mutate_epi(sim, ir = (si.flow / s.num) * 100)
# 
# # These new variables appear in the data.frame output
# as.data.frame(sim)

## ----tracker-epiby, results = "hide", message = FALSE-------------------------
set.seed(0)

nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "risk", rep(0:1, each = 50))

est <- netest(
  nw, formation = ~edges + nodefactor("risk"),
  target.stats = c(40, 20),
  coef.diss = dissolution_coefs(~offset(edges), 20),  # duration = 20
  verbose = FALSE
)

param <- param.net(inf.prob = 0.3, act.rate = 1)
init <- init.net(i.num = 5)
control <- control.net(
  type = "SI", nsims = 1, nsteps = 50,
  epi.by = "risk",
  verbose = FALSE
)

sim <- netsim(est, param, init, control)

## -----------------------------------------------------------------------------
d <- as.data.frame(sim)
names(d)

## ----tracker-example----------------------------------------------------------
epi_s_num <- function(dat) {
  needed_attributes <- c("status")
  output <- with(get_attr_list(dat, needed_attributes), {
    sum(status == "s", na.rm = TRUE)
  })
  return(output)
}

## ----tracker-commented--------------------------------------------------------
epi_prop_infected <- function(dat) {
  needed_attributes <- c("status", "active")
  output <- with(get_attr_list(dat, needed_attributes), {
    pop <- active == 1
    infected <- status == "i"
    sum(infected & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
  })
  return(output)
}

## ----tracker-list-------------------------------------------------------------
some.trackers <- list(
  prop_infected = epi_prop_infected,
  s_num         = epi_s_num
)

control <- control.net(
  type = "SI",
  nsims = 1,
  nsteps = 50,
  verbose = FALSE,
  .tracker.list = some.trackers
)

param <- param.net(
  inf.prob = 0.3,
  act.rate = 0.1
)

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),  # duration = 10, closed population
  verbose = FALSE
)

init <- init.net(i.num = 10)
sim <- netsim(est, param, init, control)

d <- as.data.frame(sim)

knitr::kable(tail(d, n = 15))

## ----record-raw, eval = FALSE-------------------------------------------------
# introspect_module <- function(dat, at) {
#   age <- get_attr(dat, "age")
# 
#   if (mean(age, na.rm = TRUE) > 50) {
#     obj <- data.frame(
#         age = age,
#         status = get_attr(dat, "status")
#     )
#     dat <- record_raw_object(dat, at, "old pop", obj)
#   }
# 
#   return(dat)
# }

## ----record-raw-access, eval = FALSE------------------------------------------
# sim <- netsim(est, param, init, control)
# 
# # Access raw records for simulation 1
# # Returns a nested list keyed by label, then time step
# raw <- sim[["raw.records"]][[1]]
# 
# # Access the "old pop" data.frame recorded at step 75
# raw[["old pop"]][[75]]

