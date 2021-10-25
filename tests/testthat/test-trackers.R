context("Custom Trackers")

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
  verbose = FALSE
)

init <- init.net(i.num = 10)
# Edges + nodematch, one-mode, closed

test_that("netsim, SI, custom trackers", {
  epi_s_num <- function(dat, at) {
    needed_attributes <- c("status")
    output <- with(get_attr_list(dat, needed_attributes), {
      out <- sum(status == "s", na.rm = TRUE)

      out
    })

    return(output)
  }

  epi_prop_infected <- function(dat, at) {
    # we need two attributes for our calculation: `status` and `active`
    needed_attributes <- c("status", "active")

    # we use `with` to simplify code
    output <- with(get_attr_list(dat, needed_attributes), {
      pop <- active == 1 # we only look at active nodes
      cond <- status == "i" # we want to know which are infected (status == "i")

      # how many are `infected` among the `active`
      out <- sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)

      out
    })

    return(output)
  }

  some.trackers <- list(
    prop_infected = epi_prop_infected,
    s_num = epi_s_num
  )

  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    trackers.FUN = trackers.net,
    tracker.list = some.trackers
  )

  param <- param.net(
    inf.prob = 0.3,
    act.rate = 0.1
  )

  mod <- netsim(est, param, init, control)

  d <- as.data.frame(mod)

  expect_true(all(c("prop_infected", "s_num") %in% names(d)))
  expect_is(d[["prop_infected"]], "numeric")
})
