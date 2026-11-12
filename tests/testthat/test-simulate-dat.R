context("simulate_dat")

# simulate_dat() calls ergm::simulate_formula() directly. These tests pin it to
# simulate() on the same formula and basis, which ends up in the same generic
# after copying the basis network on the way there.

simulate_dat_setup <- function(tergmLite, duration = 10) {
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
                         tergmLite = tergmLite, resimulate.network = TRUE,
                         verbose = FALSE)
  crosscheck.net(est, param, init, control)
  control <- netsim_validate_control(control)
  netsim_initialize(est, param, init, control, s = 1)
}

# The call simulate_dat() made before, for a layer with duration above 1
simulate_ref <- function(dat, at, time_offset, output) {
  nwparam <- get_nwparam(dat, network = 1)
  suppressWarnings(simulate(
    ~Form(nwparam$formation) + Persist(nwparam$coef.diss$dissolution),
    coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
    basis = get_network(x = dat, network = 1),
    constraints = nwparam$constraints,
    time.start = at - time_offset,
    time.offset = time_offset,
    time.slices = 1L,
    output = output,
    control = get_network_control(dat, 1, "set.control.tergm"),
    monitor = NULL,
    dynamic = TRUE
  ))
}

test_that("simulate_dat matches simulate() on a tergmLite network", {
  skip_on_cran()
  set.seed(1)
  dat <- simulate_dat_setup(tergmLite = TRUE)

  set.seed(5)
  out <- simulate_dat(dat, at = 2L)
  after <- runif(1)
  set.seed(5)
  ref <- set_network(dat, network = 1,
                     nw = simulate_ref(dat, at = 2L, time_offset = 1L, output = "final"))

  # the same draws, the same edgelist and network attributes
  expect_identical(runif(1), after)
  expect_identical(out$run$el, ref$run$el)
  expect_identical(out$run$net_attr, ref$run$net_attr)
  # and the step did something
  expect_false(identical(out$run$el[[1]], dat$run$el[[1]]))
})

test_that("simulate_dat matches simulate() on a networkDynamic", {
  skip_on_cran()
  set.seed(1)
  dat <- simulate_dat_setup(tergmLite = FALSE)

  set.seed(5)
  out <- simulate_dat(dat, at = 2L)
  after <- runif(1)
  set.seed(5)
  ref <- set_network(dat, network = 1,
                     nw = simulate_ref(dat, at = 2L, time_offset = 0L,
                                       output = "networkDynamic"))

  expect_identical(runif(1), after)
  expect_s3_class(out$run$nw[[1]], "networkDynamic")
  expect_identical(out$run$nw, ref$run$nw)
})

test_that("simulate_dat matches simulate() on a one-step layer", {
  skip_on_cran()
  set.seed(1)
  dat <- simulate_dat_setup(tergmLite = TRUE, duration = 1)
  nwparam <- get_nwparam(dat, network = 1)
  simulation_control <- get_network_control(dat, 1, "set.control.tergm")
  simulation_control$MCMC.prop.args <- list(discordance_fraction = 0)

  set.seed(5)
  out <- simulate_dat(dat, at = 2L)
  after <- runif(1)
  set.seed(5)
  ref_nw <- suppressWarnings(simulate(
    nwparam$formation,
    coef = nwparam$coef.form,
    basis = get_network(x = dat, network = 1),
    constraints = nwparam$constraints,
    time.start = 1L,
    time.offset = 1L,
    time.slices = 1L,
    output = "final",
    control = simulation_control,
    monitor = NULL,
    dynamic = TRUE
  ))
  ref <- set_network(dat, network = 1, nw = ref_nw)

  expect_identical(runif(1), after)
  expect_identical(out$run$el, ref$run$el)
})
