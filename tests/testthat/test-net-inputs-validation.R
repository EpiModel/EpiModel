context("net.inputs.R: deprecated argument guards")

# These argument names were removed in v2.6.1 (#992 / #989) and must keep
# producing hard errors so users with old scripts get a clear pointer to
# the current name. Not skipped on CRAN: the guards are cheap and the API
# contract is what they enforce.

test_that("param.net rejects the .m2 suffix", {
  expect_error(param.net(inf.prob = 0.3, inf.prob.m2 = 0.2),
               "\\.m2 suffix have been removed")
})

test_that("init.net rejects the .m2 suffix", {
  expect_error(init.net(i.num = 10, i.num.m2 = 5),
               "\\.m2 suffix have been removed")
})

test_that("control.net rejects births.FUN / deaths.FUN / depend", {
  noop <- function(dat, at) dat

  expect_error(control.net(type = "SI", nsteps = 5, births.FUN = noop),
               "births.FUN parameter has been removed")
  expect_error(control.net(type = "SI", nsteps = 5, deaths.FUN = noop),
               "deaths.FUN parameter has been removed")
  expect_error(control.net(type = "SI", nsteps = 5, depend = TRUE),
               "depend parameter has been removed")
})

context("net.inputs.R: constructor input validation")

test_that("update_params validates its inputs", {
  expect_error(update_params(list(inf.prob = 0.3), list(inf.prob = 0.5)),
               "x should be object of class param.net")
  expect_error(update_params(param.net(inf.prob = 0.3), c(inf.prob = 0.5)),
               "new.param.list should be object of class list")
})

test_that("param_random rejects mis-sized probability vectors", {
  expect_error(param_random(values = 1:3, prob = c(0.5, 0.5)),
               "incorrect number of probabilites")
})

test_that("init.net rejects conflicting initial-condition arguments", {
  expect_error(init.net(i.num = 10, status.vector = rep("s", 10)),
               "i.num OR status.vector")
  expect_error(init.net(status.vector = NULL, infTime.vector = c(1, 2, 3)),
               "infTime.vector may only be used if status.vector is used")
  expect_error(init.net(status.vector = c("s", "i", "s"),
                        infTime.vector = c(1, 2)),
               "Length of infTime.vector must match length of status.vector")
})

test_that("control.net rejects type with user modules", {
  noop <- function(dat, at) dat
  # User-supplied .FUN argument is detected as a user module; type must
  # be NULL.
  expect_error(control.net(type = "SI", nsteps = 5, custom.FUN = noop),
               "type.*null if any user")
})

test_that("control.net rejects multi-element epi.by", {
  expect_error(control.net(type = "SI", nsteps = 5,
                           epi.by = c("group", "race")),
               "epi.by currently limited to 1")
})

test_that("param.net warns on act.rate.g2", {
  expect_warning(param.net(inf.prob = 0.3, act.rate.g2 = 1),
                 "act.rate.g2.*only act.rate parameter will apply")
})

test_that("generate_random_params validates random.params structure", {
  # random.params validation runs at simulation time, inside
  # generate_random_params(). Calling that function directly lets us hit
  # the guards without spinning up a full netsim.
  p_not_list <- param.net(inf.prob = 0.3, random.params = "not a list")
  expect_error(generate_random_params(p_not_list),
               "random.params.*must be.*list")

  p_unnamed <- param.net(inf.prob = 0.3,
                         random.params = list(param_random(1:3),
                                              foo = param_random(1:3)))
  expect_error(generate_random_params(p_unnamed),
               "must be named")
})

context("crosscheck.net: end-to-end netsim input validation")

# crosscheck.net runs from inside netsim(). These are end-to-end so they
# need a netest fit -- skip on CRAN.

build_small_est <- function(two_group = FALSE) {
  nw <- network_initialize(n = 20)
  if (two_group) {
    nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 10))
  }
  netest(nw,
         formation = ~edges,
         target.stats = 5,
         coef.diss = dissolution_coefs(~offset(edges), 10, 0),
         verbose = FALSE)
}

ctrl_small <- function(type) {
  control.net(type = type, nsteps = 3, nsims = 1, verbose = FALSE)
}

test_that("netsim requires param/init/control objects of the right class", {
  skip_on_cran()
  est <- build_small_est()
  expect_error(netsim(est, list(inf.prob = 0.3), init.net(i.num = 5),
                      ctrl_small("SI")),
               "param must be an object of class param.net")
  expect_error(netsim(est, param.net(inf.prob = 0.3),
                      list(i.num = 5), ctrl_small("SI")),
               "init must be an object of class init.net")
})

test_that("crosscheck.net requires rec.rate for SIR/SIS", {
  skip_on_cran()
  est <- build_small_est()
  expect_error(netsim(est, param.net(inf.prob = 0.3),
                      init.net(i.num = 5, r.num = 0),
                      ctrl_small("SIR")),
               "Specify rec.rate in param.net")
  expect_error(netsim(est, param.net(inf.prob = 0.3),
                      init.net(i.num = 5),
                      ctrl_small("SIS")),
               "Specify rec.rate in param.net")
})

test_that("crosscheck.net requires r.num for SIR", {
  skip_on_cran()
  est <- build_small_est()
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3, rec.rate = 0.05),
                      init.net(i.num = 5),
                      ctrl_small("SIR")),
               "Specify r.num in init.net")
})

test_that("crosscheck.net validates status.vector length and values", {
  skip_on_cran()
  est <- build_small_est()
  # Length mismatch: network is size 20, status.vector is length 5.
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3),
                      init.net(status.vector = rep("s", 5)),
                      ctrl_small("SI")),
               "status.vector is unequal to size")
  # Bad values for SI (only s, i are allowed).
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3),
                      init.net(status.vector = c(rep("s", 19), "x")),
                      ctrl_small("SI")),
               "values other than .*s.* and .*i")
  # Bad values for SIR (s, i, r are allowed).
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3, rec.rate = 0.05),
                      init.net(status.vector = c(rep("s", 19), "x")),
                      ctrl_small("SIR")),
               "values other than")
})

test_that("crosscheck.net requires i.num.g2 in two-group models", {
  skip_on_cran()
  est <- build_small_est(two_group = TRUE)
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3, inf.prob.g2 = 0.2),
                      init.net(i.num = 3),
                      ctrl_small("SI")),
               "Specify i.num.g2")
})

test_that("crosscheck.net requires vital-dynamics params for two-group models", {
  skip_on_cran()
  est <- build_small_est(two_group = TRUE)
  # a.rate set without a.rate.g2 triggers the a.rate.g2 missing error.
  expect_error(netsim(est,
                      param.net(inf.prob = 0.3, inf.prob.g2 = 0.2,
                                a.rate = 0.001, ds.rate = 0.001,
                                di.rate = 0.001),
                      init.net(i.num = 3, i.num.g2 = 3),
                      ctrl_small("SI")),
               "Specify a.rate.g2")
})

test_that("crosscheck.net rejects bad x with start == 1", {
  skip_on_cran()
  expect_error(netsim("not a netest",
                      param.net(inf.prob = 0.3),
                      init.net(i.num = 5),
                      ctrl_small("SI")),
               "x must be either an object of class netest")
})

test_that("crosscheck.net validates restart inputs (start > 1)", {
  skip_on_cran()
  # start > 1 requires x to be a netsim object.
  control_restart <- control.net(type = "SI", nsteps = 5, nsims = 1,
                                 start = 3, verbose = FALSE)
  expect_error(netsim(build_small_est(),
                      param.net(inf.prob = 0.3),
                      init.net(i.num = 5),
                      control_restart),
               "x must be a netsim object")
})
