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
  expect_null(tr$control$epi.start)
  expect_equal(as.data.frame(tr, out = "mean")$time, 1:11)

  tr2 <- truncate_sim(mod, at = 10, reset.time = FALSE)
  expect_equal(tr2$control$epi.start, 10)
  df <- as.data.frame(mod)
  df2 <- as.data.frame(tr2)
  expect_equal(df2$time, rep(10:20, 2))
  expect_equal(df2$i.num, df$i.num[df$time >= 10])
})

test_that("truncate_sim.icm reads `at` on the time steps of a truncated object", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1)
  init <- init.icm(s.num = 50, i.num = 2)
  control <- control.icm(type = "SIS", nsteps = 20, nsims = 2,
                         verbose = FALSE)
  mod <- icm(param, init, control)
  tr <- truncate_sim(mod, at = 10, reset.time = FALSE)

  expect_equal(truncate_sim(tr, at = 15, reset.time = FALSE),
               truncate_sim(mod, at = 15, reset.time = FALSE))
  expect_equal(truncate_sim(tr, at = 15), truncate_sim(mod, at = 15))

  expect_error(truncate_sim(tr, at = 9), "between 10 and 20")
  expect_error(truncate_sim(tr, at = 21), "between 10 and 20")
  expect_error(truncate_sim(mod, at = 0), "between 1 and 20")
  expect_error(truncate_sim(mod, at = c(5, 6)), "between 1 and 20")
})

test_that("plot, summary and comp_plot read the time steps of an icm truncated without reset", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1)
  init <- init.icm(s.num = 50, i.num = 2)
  control <- control.icm(type = "SIS", nsteps = 20, nsims = 2,
                         verbose = FALSE)
  mod <- icm(param, init, control)
  tr <- truncate_sim(mod, at = 10, reset.time = FALSE)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(tr, y = "i.num")
  expect_equal(par("usr")[1:2], extendrange(c(9, 20), f = 0.04))

  expect_output(summary(tr, at = 15), "Time: 15")
  expect_error(summary(tr, at = 5), "Specify a time step between 10 and 20")
  expect_silent(comp_plot(tr, at = 15))
  expect_error(comp_plot(tr, at = 5), "Specify a timestep between 10 and 20")
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

build_restart_sim <- function(nsteps = 5, nsims = 1, inf.prob = 0.3) {
  nw <- network_initialize(n = 30)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  # Set resimulate.network = TRUE explicitly so control.net() doesn't emit
  # the "tergmLite = TRUE ... resetting resimulate.network" notice.
  control <- control.net(type = "SI", nsteps = nsteps, nsims = nsims,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  netsim(est, param.net(inf.prob = inf.prob),
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

test_that("make_restart_point accepts a simulation without transmissions", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 6, inf.prob = 0)
  expect_silent(rp <- make_restart_point(mod, time_attrs = c(), keep_steps = 1))
  expect_equal(nrow(get_transmat(rp)), 0)
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
  # seed with transmissions before and within the kept steps of `sim2`
  set.seed(4)
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

  # `nwstats` and `transmat` are trimmed
  expect_equal(
    as.data.frame(rp$stats$nwstats$sim1[[1]]),
    as.data.frame(tail(mod$stats$nwstats$sim2[[1]], 2)),
    ignore_attr = TRUE
  )
  tm <- mod$stats$transmat$sim2
  kept <- tm$at > 6 - 2
  expect_true(any(kept) && !all(kept))
  expect_equal(rp$stats$transmat$sim1$at, tm$at[kept] - (6 - 2))
})

test_that("make_restart_point re-bases the unique IDs of the transmission matrix", {
  skip_on_cran()
  set.seed(12)
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.05),
                verbose = FALSE)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  # high departure rates so that the lowest unique IDs depart
  param <- param.net(inf.prob = 0.5, act.rate = 2,
                     a.rate = 0.3, ds.rate = 0.3, di.rate = 0.3)
  mod <- netsim(est, param, init.net(i.num = 10), control)

  uid_offset <- min(mod$run$sim1$attr$unique_id) - 1
  expect_gt(uid_offset, 0)

  rp <- make_restart_point(mod, time_attrs = "infTime", keep_steps = 10)
  tm <- mod$stats$transmat$sim1
  expect_equal(rp$stats$transmat$sim1$sus, tm$sus - uid_offset)
  expect_equal(rp$stats$transmat$sim1$inf, tm$inf - uid_offset)
})

test_that("make_restart_point keeps the empty transmission matrix of a simulation without transmissions", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 4)
  # what `netsim` saves for a simulation without transmissions
  mod$stats$transmat$sim1 <- dplyr::tibble()
  expect_silent(
    rp <- make_restart_point(mod, time_attrs = "infTime", keep_steps = 2)
  )
  expect_equal(rp$stats$transmat$sim1, dplyr::tibble())
})

test_that("simulations restarted from make_restart_point keep nwstats and transmat in line with epi", {
  skip_on_cran()
  # seed with transmissions before and within the kept steps
  set.seed(4)
  mod <- build_restart_sim(nsteps = 6)
  rp <- make_restart_point(mod, time_attrs = "infTime", keep_steps = 3)
  rp_tm <- rp$stats$transmat$sim1
  expect_gt(nrow(rp_tm), 0)

  control <- control.net(type = "SI", start = 4, nsteps = 8, nsims = 2,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  mod2 <- netsim(rp, param.net(inf.prob = 0.3), init.net(i.num = 5), control)

  # the epi trackers, the network statistics and the transmission matrix hold
  # the steps kept in the restart point and the new ones
  expect_equal(nrow(mod2$epi$i.num), 8)
  for (s in 1:2) {
    nwstats2 <- mod2$stats$nwstats[[s]][[1]]
    expect_equal(nrow(nwstats2), 8)
    expect_false(anyNA(nwstats2))
    expect_equal(
      as.matrix(nwstats2[1:3, , drop = FALSE]),
      as.matrix(rp$stats$nwstats$sim1[[1]])
    )

    tm2 <- mod2$stats$transmat[[s]]
    expect_equal(
      unname(as.matrix(tm2[seq_len(nrow(rp_tm)), names(rp_tm), drop = FALSE])),
      unname(as.matrix(rp_tm))
    )
    expect_true(all(tm2$at[-seq_len(nrow(rp_tm))] >= 4))
  }

  nws <- get_nwstats(mod2)
  expect_equal(nws$time, rep(1:8, 2))
  expect_equal(nws$sim, rep(1:2, each = 8))
  expect_equal(as.data.frame(mod2)$time, nws$time)

  rp2 <- make_restart_point(mod2, time_attrs = "infTime", keep_steps = 2)
  expect_equal(
    as.data.frame(rp2$stats$nwstats$sim1[[1]]),
    as.data.frame(tail(mod2$stats$nwstats$sim1[[1]], 2)),
    ignore_attr = TRUE
  )
  tm2 <- mod2$stats$transmat$sim1
  expect_equal(rp2$stats$transmat$sim1$at, tm2$at[tm2$at > 8 - 2] - (8 - 2))
})

test_that("restarted simulations label their time steps from 1", {
  skip_on_cran()
  set.seed(1)
  mod <- build_restart_sim(nsteps = 5)
  control <- control.net(type = "SI", start = 6, nsteps = 10, nsims = 2,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)
  mod2 <- netsim(mod, param.net(inf.prob = 0.3), init.net(i.num = 5), control)

  # the epi trackers hold the steps of `x` followed by the new ones, in line
  # with the network statistics
  df <- as.data.frame(mod2)
  expect_equal(df$time, rep(1:10, 2))
  expect_equal(df$i.num[df$sim == 2 & df$time <= 5], mod$epi$i.num$sim1)
  expect_equal(get_nwstats(mod2)$time, df$time)
  df_mean <- as.data.frame(mod2, out = "mean")
  expect_equal(df_mean$time, 1:10)

  # truncating at a new step without reset keeps its label
  tr <- truncate_sim(mod2, at = 6, reset.time = FALSE)
  df_tr <- as.data.frame(tr)
  expect_equal(df_tr$time, rep(6:10, 2))
  expect_equal(df_tr$i.num, df$i.num[df$time >= 6])

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(mod2, y = "i.num")
  expect_equal(par("usr")[1:2], extendrange(c(0, 10), f = 0.04))
  plot(tr, y = "i.num")
  expect_equal(par("usr")[1:2], extendrange(c(5, 10), f = 0.04))
  plot(df_tr, y = "i.num")
  expect_equal(par("usr")[1:2], extendrange(c(5, 10), f = 0.04))

  # summary and comp_plot report the step asked for
  summary_i_num <- function(x, at) {
    out <- capture.output(summary(x, at = at))
    as.numeric(strsplit(grep("^Infect", out, value = TRUE), " +")[[1]][2])
  }
  expect_equal(summary_i_num(mod2, at = 8), round(df_mean$i.num[8], 3))
  expect_equal(summary_i_num(tr, at = 8), round(df_mean$i.num[8], 3))
  expect_error(summary(tr, at = 5), "Specify at between 6 and 10")
  expect_silent(comp_plot(mod2, at = 8))
})

test_that("restarting from a simulation without transmissions keeps the transmission matrix a data frame", {
  skip_on_cran()
  set.seed(2)
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  init <- init.net(i.num = 25)
  control <- control.net(type = "SI", nsteps = 3, tergmLite = TRUE,
                         resimulate.network = TRUE, save.run = TRUE,
                         verbose = FALSE)
  mod <- netsim(est, param.net(inf.prob = 0), init, control)
  expect_equal(mod$stats$transmat$sim1, dplyr::tibble())

  control$start <- 4
  control$nsteps <- 6

  # no transmissions after the restart either
  expect_silent(mod2 <- netsim(mod, param.net(inf.prob = 0), init, control))
  expect_equal(mod2$stats$transmat$sim1, dplyr::tibble())

  # transmissions after the restart
  expect_silent(mod3 <- netsim(mod, param.net(inf.prob = 1), init, control))
  tm <- mod3$stats$transmat$sim1
  expect_gt(nrow(tm), 0)
  expect_false(inherits(tm, "tbl_df"))
  expect_true(all(tm$at >= 4))
})

test_that("restarting from x without nwstats records the earlier steps as NA", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, tergmLite = TRUE,
                         resimulate.network = TRUE, save.run = TRUE,
                         save.nwstats = FALSE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_null(x$stats$nwstats)

  rp <- make_restart_point(x, time_attrs = "infTime")
  for (x_rs in list(x, rp)) {
    start <- x_rs$control$nsteps + 1
    nsteps <- start + 4
    for (nsims in 1:2) {
      control_rs <- control.net(type = "SI", start = start, nsteps = nsteps,
                                nsims = nsims, tergmLite = TRUE,
                                resimulate.network = TRUE, save.run = TRUE,
                                verbose = FALSE)
      mod <- netsim(x_rs, param, init, control_rs)

      expect_equal(nrow(mod$epi$i.num), nsteps)
      for (s in seq_len(nsims)) {
        nwstats <- mod$stats$nwstats[[s]][[1]]
        expect_equal(nrow(nwstats), nsteps)
        expect_true(all(is.na(nwstats[seq_len(start - 1), ])))
        expect_false(anyNA(nwstats[start:nsteps, ]))
      }

      nws <- get_nwstats(mod)
      expect_equal(nws$time, rep(seq_len(nsteps), nsims))
      expect_equal(nws$sim, rep(seq_len(nsims), each = nsteps))
    }
  }
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

test_that("netsim refuses to restart from objects of another EpiModel version", {
  skip_on_cran()
  mod <- build_restart_sim(nsteps = 3, nsims = 2)
  pkg_ver <- utils::packageVersion("EpiModel")
  expect_equal(attr(mod, "version"), pkg_ver)

  # the version is kept by the functions preparing a restart
  rp <- make_restart_point(mod, time_attrs = "infTime", sim_num = 2)
  expect_equal(attr(get_sims(mod, sims = 1), "version"), pkg_ver)
  expect_equal(attr(merge(mod, mod), "version"), pkg_ver)
  expect_equal(attr(rp, "version"), pkg_ver)

  control <- control.net(type = "SI", start = 2, nsteps = 4, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  restart <- function(x) {
    netsim(x, param.net(inf.prob = 0.3), init.net(i.num = 5), control)
  }
  expect_s3_class(restart(rp), "netsim")

  # objects from 2.6.x or earlier have no version
  attr(rp, "version") <- NULL
  expect_error(restart(rp), "comes from version v(<2.7)", fixed = TRUE)

  attr(rp, "version") <- package_version("2.6.2")
  expect_error(restart(rp), "comes from version v2.6.2", fixed = TRUE)
})
