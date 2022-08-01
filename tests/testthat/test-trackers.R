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
      sum(status == "s", na.rm = TRUE)
    })
    return(output)
  }

  epi_prop_infected <- function(dat, at) {
    needed_attributes <- c("status", "active")
    output <- with(get_attr_list(dat, needed_attributes), {
      pop <- active == 1
      cond <- status == "i"
      sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
    })
    return(output)
  }

  some.trackers <- list(
    prop_infected = epi_prop_infected,
    s_num = epi_s_num
  )

  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    .tracker.list = some.trackers
  )

  param <- param.net(
    inf.prob = 0.3,
    act.rate = 0.1
  )

  mod <- netsim(est, param, init, control)

  d <- as.data.frame(mod)

  expect_true(all(c("prop_infected", "s_num") %in% names(d)))
  expect_is(d[["prop_infected"]], "numeric")
  # the custom epi trackers are not run during intialization so the first value
  # is always NA
  expect_true(all(d$s_num[2:50] == d$s.num[2:50]))
})
