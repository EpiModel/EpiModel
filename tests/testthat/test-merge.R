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
