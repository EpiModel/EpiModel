context("ergm term version")

# simulate_dat() hands ergm's term initializers the installed ergm version
# through `term.options`, which is the value they would otherwise read from
# disk as their default.

test_that("set_term_version fills in the installed ergm version", {
  control <- tergm::control.simulate.formula.tergm()
  expect_null(control$term.options$version)
  out <- set_term_version(control)
  expect_identical(out$term.options$version, utils::packageVersion("ergm"))
  expect_identical(class(out), class(control))
  expect_identical(out[names(out) != "term.options"],
                   control[names(control) != "term.options"])
})

test_that("set_term_version keeps what the user set", {
  old <- as.package_version("3.9.4")
  control <- tergm::control.simulate.formula.tergm(
    term.options = list(version = old, interact.dependent = "silent")
  )
  expect_identical(set_term_version(control), control)

  control <- tergm::control.simulate.formula.tergm(
    term.options = list(interact.dependent = "silent")
  )
  out <- set_term_version(control)
  expect_identical(out$term.options$interact.dependent, "silent")
  expect_identical(out$term.options$version, utils::packageVersion("ergm"))
})

test_that("simulate_dat gives the same network as the default version", {
  skip_on_cran()
  set.seed(1)
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "grp", rep(1:2, each = 50))
  est <- netest(
    nw,
    formation = ~edges + nodematch("grp") + nodefactor("grp"),
    target.stats = c(30, 20, 28),
    coef.diss = dissolution_coefs(~offset(edges), 10),
    verbose = FALSE
  )
  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)
  crosscheck.net(est, param, init, control)
  control <- netsim_validate_control(control)
  dat <- netsim_initialize(est, param, init, control, s = 1)
  nwparam <- get_nwparam(dat, network = 1)

  set.seed(5)
  out <- simulate_dat(dat, at = 2L)
  after <- runif(1)
  # the same call with the control as the user gave it, no version in it
  set.seed(5)
  ref_nw <- suppressWarnings(simulate(
    ~Form(nwparam$formation) + Persist(nwparam$coef.diss$dissolution),
    coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
    basis = get_network(x = dat, network = 1),
    constraints = nwparam$constraints,
    time.start = 1L,
    time.offset = 1L,
    time.slices = 1L,
    output = "final",
    control = get_network_control(dat, 1, "set.control.tergm"),
    monitor = NULL,
    dynamic = TRUE
  ))
  ref <- set_network(dat, network = 1, nw = ref_nw)

  expect_null(get_network_control(dat, 1, "set.control.tergm")$term.options$version)
  expect_identical(runif(1), after)
  expect_identical(out$run$el, ref$run$el)
  expect_false(identical(out$run$el[[1]], dat$run$el[[1]]))
})
