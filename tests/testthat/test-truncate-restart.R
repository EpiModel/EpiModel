context("truncate_sim S3 methods")

test_that("truncate_sim.dcm trims epi and resets timesteps", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05)
  init <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 20)
  mod <- dcm(param, init, control)
  full_rows <- nrow(mod$epi[[1]])

  # reset.time = TRUE (default): timesteps start back at 1.
  tr <- truncate_sim(mod, at = 10)
  expect_equal(nrow(tr$epi[[1]]), full_rows - 10 + 1)
  expect_equal(min(tr$control$timesteps), 1)
  expect_equal(tr$control$nsteps, max(tr$control$timesteps))

  # reset.time = FALSE: timesteps keep their original numbering.
  tr2 <- truncate_sim(mod, at = 10, reset.time = FALSE)
  expect_equal(nrow(tr2$epi[[1]]), full_rows - 10 + 1)
  expect_equal(min(tr2$control$timesteps), 10)
})

test_that("truncate_sim.dcm errors when `at` is not in the timesteps vector", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05)
  init <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 20)
  mod <- dcm(param, init, control)

  expect_error(truncate_sim(mod, at = 999),
               "at is not in the control\\$timesteps vector")
})

test_that("truncate_sim.icm trims epi and resets start", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1)
  init <- init.icm(s.num = 50, i.num = 2)
  control <- control.icm(type = "SIS", nsteps = 20, nsims = 2,
                         verbose = FALSE)
  mod <- icm(param, init, control)

  tr <- truncate_sim(mod, at = 10)
  expect_equal(nrow(tr$epi[[1]]), 20 - 10 + 1)
  expect_equal(tr$control$start, 1)

  tr2 <- truncate_sim(mod, at = 10, reset.time = FALSE)
  expect_equal(tr2$control$start, 10)
})

test_that("truncate_sim.netsim delegates to the icm method", {
  skip_on_cran()
  nw <- network_initialize(n = 30)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  control <- control.net(type = "SI", nsteps = 12, nsims = 1,
                         verbose = FALSE)
  mod <- netsim(est, param.net(inf.prob = 0.3),
                init.net(i.num = 5), control)

  tr <- truncate_sim(mod, at = 6)
  expect_equal(nrow(tr$epi[[1]]), 12 - 6 + 1)
  expect_equal(tr$control$start, 1)
})

context("make_restart_point validation and roundtrip")

build_restart_sim <- function(nsteps = 5, nsims = 1) {
  nw <- network_initialize(n = 30)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  # Set resimulate.network = TRUE explicitly so control.net() doesn't emit
  # the "tergmLite = TRUE ... resetting resimulate.network" notice.
  control <- control.net(type = "SI", nsteps = nsteps, nsims = nsims,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  netsim(est, param.net(inf.prob = 0.3),
         init.net(i.num = 5), control)
}

test_that("make_restart_point rejects non-netsim input", {
  expect_error(make_restart_point(list(a = 1), time_attrs = c()),
               "must be.*netsim")
})

test_that("make_restart_point rejects non-tergmLite simulations", {
  skip_on_cran()
  nw <- network_initialize(n = 30)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  # tergmLite = FALSE: should not be accepted.
  control <- control.net(type = "SI", nsteps = 5, nsims = 1,
                         tergmLite = FALSE, save.run = TRUE,
                         verbose = FALSE)
  mod <- netsim(est, param.net(inf.prob = 0.3),
                init.net(i.num = 5), control)

  expect_error(make_restart_point(mod, time_attrs = c()),
               "tergmLite == TRUE")
})

test_that("make_restart_point validates sim_num and keep_steps ranges", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 5)

  expect_error(make_restart_point(mod, time_attrs = c(), sim_num = 0),
               "sim_num.*>= 1.*<=")
  expect_error(make_restart_point(mod, time_attrs = c(), sim_num = 99),
               "sim_num.*>= 1.*<=")
  expect_error(make_restart_point(mod, time_attrs = c(), sim_num = c(1, 2)),
               "sim_num.*single")
  expect_error(make_restart_point(mod, time_attrs = c(), keep_steps = 0),
               "keep_steps.*>= 1.*<=")
  expect_error(make_restart_point(mod, time_attrs = c(), keep_steps = 99),
               "keep_steps.*>= 1.*<=")
})

test_that("make_restart_point produces a trimmed netsim with reset time", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 6)

  rp <- make_restart_point(mod, time_attrs = c(), keep_steps = 1)

  expect_s3_class(rp, "netsim")
  expect_equal(rp$control$start, 1)
  expect_equal(rp$control$nsteps, 1)
  expect_equal(nrow(rp$epi[[1]]), 1)
  # unique_id offset: smallest id should be 1 after trimming.
  expect_equal(min(rp$run$sim1$attr$unique_id), 1)
  # Default-trimmed extras.
  expect_equal(rp$attr.history, list())
  expect_equal(rp$raw.records, list())
})

test_that("make_restart_point keeps multiple history steps when requested", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 6)

  rp <- make_restart_point(mod, time_attrs = c(), keep_steps = 3)
  expect_equal(rp$control$nsteps, 3)
  expect_equal(nrow(rp$epi[[1]]), 3)
})

test_that("make_restart_point reports missing time attributes by name", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 5)

  expect_error(make_restart_point(mod, time_attrs = "not.an.attr"),
               "not.an.attr")
})

test_that("make_restart_point extracts and trims the requested simulation", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 6, nsims = 3)

  rp <- make_restart_point(mod, time_attrs = c(), sim_num = 2, keep_steps = 2)

  # The restart point always holds a single simulation, renamed `sim1`
  expect_equal(rp$control$nsims, 1)
  expect_named(rp$run, "sim1")
  expect_named(rp$coef.form, "sim1")

  # The kept simulation is the requested one
  expect_equal(rp$run$sim1$attr$status, mod$run$sim2$attr$status)
  expect_equal(rp$coef.form$sim1, mod$coef.form$sim2)
  expect_equal(rp$epi$i.num[, 1], tail(mod$epi$i.num[, 2], 2))

  # UIDs and time are re-based
  expect_equal(min(rp$run$sim1$attr$unique_id), 1)
  expect_equal(
    rp$run$sim1$attr$entrTime,
    mod$run$sim2$attr$entrTime - (6 - 2)
  )

  # `nwstats` and `transmat` are not kept
  expect_false(is.null(mod$stats$nwstats))
  expect_false(is.null(mod$stats$transmat))
  expect_null(rp$stats$nwstats)
  expect_null(rp$stats$transmat)
})

test_that("simulations restarted from make_restart_point record only their new steps", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 6)
  rp <- make_restart_point(mod, time_attrs = "infTime", keep_steps = 3)

  control <- control.net(type = "SI", start = 4, nsteps = 8,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  mod2 <- netsim(rp, param.net(inf.prob = 0.3), init.net(i.num = 5), control)

  # the epi trackers hold the steps kept in the restart point and the new ones,
  # the network statistics and the transmission matrix only the new ones
  expect_equal(nrow(mod2$epi$i.num), 8)
  nwstats2 <- mod2$stats$nwstats$sim1[[1]]
  expect_equal(nrow(nwstats2), 8 - 4 + 1)
  expect_false(anyNA(nwstats2))
  expect_true(all(mod2$stats$transmat$sim1$at >= 4))

  rp2 <- make_restart_point(mod2, time_attrs = "infTime", keep_steps = 2)
  expect_null(rp2$stats$nwstats)
  expect_null(rp2$stats$transmat)
})

test_that("restarted simulations all restart from the first simulation of x", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 5, nsims = 2)
  # tag each simulation of `x` in an extra epi column
  mod$epi$tag <- data.frame(sim1 = rep(1, 5), sim2 = rep(2, 5))

  control <- control.net(type = "SI", start = 6, nsteps = 7, nsims = 3,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  expect_warning(
    mod2 <- netsim(mod, param.net(inf.prob = 0.3), init.net(i.num = 5), control),
    "holds 2 simulations: every simulation restarts from the first one"
  )

  expect_equal(unname(unlist(mod2$epi$tag[1, ])), c(1, 1, 1))
  for (s in 1:3) {
    expect_equal(mod2$epi$i.num[1:5, s], mod$epi$i.num[1:5, 1])
  }

  # a single simulation selected with `get_sims()` restarts without warning
  expect_silent(
    mod3 <- netsim(get_sims(mod, sims = 2), param.net(inf.prob = 0.3),
                   init.net(i.num = 5), control)
  )
  expect_equal(unname(unlist(mod3$epi$tag[1, ])), c(2, 2, 2))
})

test_that("restarted simulations restore the groups parameter from run", {
  skip_on_cran()
  nw <- network_initialize(n = 30)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 15))
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.3)
  init <- init.net(i.num = 3, i.num.g2 = 3)
  control <- control.net(type = "SI", nsteps = 5, tergmLite = TRUE,
                         resimulate.network = TRUE, save.run = TRUE,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_equal(mod$run$sim1$groups, 2)

  # `groups` comes from `run`, as neither `x` nor `param` provide it
  rp <- make_restart_point(mod, time_attrs = "infTime")
  rp$param <- NULL
  control$start <- 2
  control$nsteps <- 4
  mod2 <- netsim(rp, param, init, control)
  expect_equal(mod2$param$groups, 2)
  expect_false(anyNA(mod2$epi$num.g2))
})
