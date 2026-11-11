context("simulate_dat edgelist")

# simulate_dat() reads a tergmLite layer's edgelist from the sampler's state.
# These tests pin it to the earlier route through tergm's final network.

state_edgelist_setup <- function(duration = 10, track.duration = FALSE) {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "grp", rep(1:2, each = 50))
  est <- netest(
    nw,
    formation = ~edges + nodematch("grp"),
    target.stats = c(30, 20),
    coef.diss = dissolution_coefs(~offset(edges), duration),
    verbose = FALSE
  )
  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         tergmLite.track.duration = track.duration,
                         verbose = FALSE)
  crosscheck.net(est, param, init, control)
  control <- netsim_validate_control(control)
  netsim_initialize(est, param, init, control, s = 1)
}

# The final network of the step, as simulate_dat() asked tergm for it before
simulate_final <- function(dat, at) {
  nwparam <- get_nwparam(dat, network = 1)
  simulation_control <- get_network_control(dat, 1, "set.control.tergm")
  if (nwparam$coef.diss$duration[1] > 1) {
    formula <- ~Form(nwparam$formation) +
      Persist(nwparam$coef.diss$dissolution)
    coef <- c(nwparam$coef.form, nwparam$coef.diss$coef.adj)
  } else {
    formula <- nwparam$formation
    coef <- nwparam$coef.form
    simulation_control$MCMC.prop.args <- list(discordance_fraction = 0)
  }
  suppressWarnings(simulate(
    formula,
    coef = coef,
    basis = get_network(x = dat, network = 1),
    constraints = nwparam$constraints,
    time.start = at - 1L,
    time.offset = 1L,
    time.slices = 1L,
    output = "final",
    control = simulation_control,
    monitor = NULL,
    dynamic = TRUE
  ))
}

expect_same_step <- function(dat) {
  set.seed(5)
  out <- simulate_dat(dat, at = 2L)
  after <- runif(1)
  set.seed(5)
  ref <- set_network(dat, network = 1, nw = simulate_final(dat, at = 2L))

  # the same draws, the same edgelist and network attributes
  expect_identical(runif(1), after)
  expect_identical(out$run$el, ref$run$el)
  expect_identical(out$run$net_attr, ref$run$net_attr)
  # and the step did something
  expect_false(identical(out$run$el[[1]], dat$run$el[[1]]))
  invisible(out)
}

test_that("simulate_dat reads the final network's edgelist from the state", {
  skip_on_cran()
  set.seed(1)
  dat <- state_edgelist_setup()
  expect_same_step(dat)
})

test_that("simulate_dat does the same on a one-step layer", {
  skip_on_cran()
  set.seed(1)
  dat <- state_edgelist_setup(duration = 1)
  expect_same_step(dat)
})

test_that("a layer that tracks durations keeps its time and lasttoggle", {
  skip_on_cran()
  set.seed(1)
  dat <- state_edgelist_setup(track.duration = TRUE)
  out <- expect_same_step(dat)
  expect_identical(out$run$net_attr[[1]][["time"]], 2L)
  expect_true(nrow(out$run$net_attr[[1]][["lasttoggle"]]) > 0)
})

test_that("ergm_state_edgelist is the edgelist of the state's network", {
  set.seed(1)
  n <- 50L
  el <- cbind(sample(1:20, 30, replace = TRUE), sample(21:50, 30, replace = TRUE))
  el <- unique(el)
  el <- el[order(el[, 1], el[, 2]), ]
  attr(el, "n") <- n
  nw <- networkLite(el, list(grp = rep(1:2, length.out = n)))
  state <- ergm::ergm_state(nw)
  expect_identical(ergm_state_edgelist(state), as.edgelist(nw))
  expect_identical(ergm_state_edgelist(state), as.edgelist(as.network(state)))

  # the state's edges in any order
  shuffled <- state
  shuffled$el <- shuffled$el[sample(nrow(shuffled$el)), ]
  expect_identical(ergm_state_edgelist(shuffled), as.edgelist(nw))

  # no edges
  empty <- el[0, , drop = FALSE]
  attr(empty, "n") <- n
  nw0 <- networkLite(empty, list(grp = rep(1:2, length.out = n)))
  state0 <- ergm::ergm_state(nw0)
  expect_identical(ergm_state_edgelist(state0), as.edgelist(nw0))
  expect_identical(nrow(ergm_state_edgelist(state0)), 0L)
})
