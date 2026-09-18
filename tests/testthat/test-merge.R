context("Stochastic Model Merging")


# merge.icm ---------------------------------------------------------------

test_that("merge for ICM", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 3, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
  expect_true(z$control$nsims == 4)
  expect_true(dim(z$epi$i.num)[2] == 4)

  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
})

test_that("merge 1 sim each", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
})

test_that("merge errors", {
  param <- param.icm(inf.prob = 0.4, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 3, verbose = FALSE)
  x <- icm(param, init, control)
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  expect_error(merge(x, y), "x and y have different parameters")
})

# merge.netsim ------------------------------------------------------------

test_that("merge for netsim", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 20, nsims = 2,
                         save.nwstats = TRUE,
                         nwstats.formula = ~edges + degree(0),
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)
  expect_is(z, "netsim")
  expect_true(z$control$nsims == 4)
  expect_true(dim(z$epi$i.num)[2] == 4)
})

test_that("merge for netsim", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 20, nsims = 2,
                         save.nwstats = TRUE,
                         nwstats.formula = ~edges + degree(0),
                         verbose = FALSE, save.other = "run")
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y, keep.other = TRUE)
  expect_is(z, "netsim")
  expect_true(length(z$run) == 4)
  expect_true(length(z$run[[1]]$attr) == 6)
  z <- merge(x, y, keep.other = FALSE)
  expect_true(any(names(z) == "run") == FALSE)
})

test_that("merge works for open sims saving nw stats", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.9, a.rate = 0.01, ds.rate = 0.01,
                     di.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 5, save.nwstats = TRUE,
                         nwstats.formula =
                           ~edges + meandeg + degree(0) + concurrent,
                         resimulate.network = TRUE, verbose = FALSE)

  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)

  nws <- get_nwstats(z)
  expect_true(nrow(nws) == 10)
  expect_true(length(unique(nws$sim)) == 2)

})

test_that("merge.netsim works as expected for transmat", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  mod <- netsim(est, param, init, control)

  expect_equal(length(mod$stats$transmat), 2)

  mod2 <- merge(mod, mod)
  expect_equal(length(mod2$stats$transmat), 4)

  mod3 <- merge(mod, mod, keep.transmat = FALSE)
  expect_true(is.null(mod3$stats$transmat))

  mod4 <- merge(mod2, mod3)
  expect_true(is.null(mod4$stats$transmat))
})

test_that("merge.netsim preserves random parameter values", {
  random.params <- list(
    act.rate = function() NULL,
    dummy.strat.param = function() NULL
  )

  make_mod <- function(act.rate, dummy.strat.param) {
    nsims <- length(act.rate)
    simnames <- paste0("sim", seq_len(nsims))
    structure(
      list(
        param = list(
          inf.prob = 0.3,
          act.rate = act.rate[1],
          dummy.strat.param = dummy.strat.param[[1]],
          random.params = random.params,
          random.params.values = list(
            act.rate = act.rate,
            dummy.strat.param = dummy.strat.param
          )
        ),
        control = list(
          nsims = nsims,
          save.other = character(0),
          monitors = NULL,
          nwstats.formula = NULL
        ),
        epi = list(
          i.num = data.frame(
            matrix(seq_len(nsims), nrow = 1, dimnames = list(NULL, simnames))
          )
        ),
        stats = list(nwstats = NULL, transmat = NULL),
        run = NULL,
        network = NULL,
        diss.stats = NULL
      ),
      class = "netsim"
    )
  }

  x <- make_mod(c(0.1, 0.2), list(c(1, 2), c(3, 4)))
  y <- make_mod(c(0.3, 0.4), list(c(5, 6), c(7, 8)))

  z <- merge(x, y)
  d.set <- get_param_set(z)

  expect_equal(z$control$nsims, 4)
  expect_equal(z$param$random.params.values$act.rate, c(0.1, 0.2, 0.3, 0.4))
  expect_equal(z$param$random.params.values$dummy.strat.param,
               list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)))
  expect_equal(d.set$act.rate, c(0.1, 0.2, 0.3, 0.4))
  expect_equal(d.set$dummy.strat.param_1, c(1, 3, 5, 7))
  expect_equal(d.set$dummy.strat.param_2, c(2, 4, 6, 8))
})

test_that("merge and print work as expected for save.other", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.other = c("run", "el"))
  mod <- netsim(est, param, init, control)

  capture_output(
    print(mod)
  )
  expect_output(print(mod), "Other Elements: run el")
  expect_equal(length(mod[["run"]]), 2)
  expect_equal(length(mod[["el"]]), 2)

  mod2 <- merge(mod, mod)
  expect_output(print(mod2), "Other Elements: run el")
  expect_equal(length(mod2[["run"]]), 4)
  expect_equal(length(mod2[["el"]]), 4)

  mod3 <- merge(mod, mod, keep.other = FALSE)
  expect_error(expect_output(print(mod3), "Other Elements"))
  expect_true(is.null(mod3[["run"]]))
  expect_true(is.null(mod3[["el"]]))
})

test_that("merge.netsim merges run, cumulative edgelist, and coef.form", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         cumulative.edgelist = TRUE,
                         save.cumulative.edgelist = TRUE,
                         save.run = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)

  simnames <- paste0("sim", 1:4)

  z <- merge(x, y, keep.cumulative.edgelist = TRUE)
  expect_equal(z$control$nsims, 4)
  expect_named(z$run, simnames)
  expect_named(z$cumulative.edgelist, simnames)
  expect_named(z$coef.form, simnames)
  expect_equal(z$run[["sim3"]], y$run[["sim1"]])
  expect_equal(z$cumulative.edgelist[["sim4"]], y$cumulative.edgelist[["sim2"]])
  expect_equal(z$coef.form[["sim3"]], y$coef.form[["sim1"]])

  # cumulative edgelists are dropped by default
  expect_null(merge(x, y)$cumulative.edgelist)

  z2 <- merge(x, y, keep.run = FALSE, keep.cumulative.edgelist = FALSE)
  expect_equal(z2$control$nsims, 4)
  expect_null(z2$run)
  expect_null(z2$cumulative.edgelist)
  # `coef.form` is always kept: it is small and needed by other accessors
  expect_named(z2$coef.form, simnames)
})

test_that("merge.netsim merges the recorded histories", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  # record one value per time step, plus a raw object
  test_logger <- function(dat, at) {
    dat <- record_attr_history(dat, at, "test.attr", get_posit_ids(dat), at)
    dat <- record_raw_object(dat, at, "test.obj", at)
    return(dat)
  }

  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 5)
  control <- control.net(type = NULL, nsteps = 5, nsims = 2, verbose = FALSE,
                         infection.FUN = infection.net,
                         logger.FUN = test_logger)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)

  simnames <- paste0("sim", 1:4)

  z <- merge(x, y)
  expect_named(z$attr.history, simnames)
  expect_named(z$raw.records, simnames)
  expect_length(z$raw.records, 4)

  # `get_attr_history` reads the simulation number off the element names, so
  # the runs coming from `y` must be reported as sims 3 and 4
  hist <- get_attr_history(z)
  expect_equal(sort(unique(hist$test.attr$sim)), c(1, 2, 3, 4))

  z2 <- merge(x, y, keep.attr.history = FALSE)
  expect_null(z2$attr.history)
  expect_null(z2$raw.records)
})

test_that("merge.netsim output can be used as a restart pool", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.run = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)

  z <- merge(x, y)

  control.rs <- control.net(type = "SI", start = 6, nsteps = 8, nsims = 4,
                            tergmLite = TRUE, resimulate.network = TRUE,
                            save.run = TRUE, verbose = FALSE)
  rs <- netsim(z, param, init, control.rs)

  expect_is(rs, "netsim")
  expect_equal(rs$control$nsims, 4)
  # each output simulation restarts from its own source run in the pool
  restart_src <- unlist(lapply(rs$run, function(r) r[["_restart_simnum"]]))
  expect_equal(unname(restart_src), 1:4)
})

test_that("merge.netsim lets save.other drive run when keep.other is FALSE", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         save.other = "run", verbose = FALSE)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)

  expect_length(merge(x, y, keep.other = TRUE)$run, 4)
  expect_null(merge(x, y, keep.other = FALSE)$run)
  # `keep.run = FALSE` does not corrupt a `run` managed by `save.other`
  expect_length(merge(x, y, keep.other = TRUE, keep.run = FALSE)$run, 4)
})
