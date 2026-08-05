context("plot sims selection for stochastic models")

# `plot.icm` and `plot_netsim_epi` are sibling paths over the same summary
# statistics, so the `sims` behavior is checked identically for both.

# Capture what a plot call draws, by mocking the two drawing functions used for
# the mean line, the simulation lines, and the quantile band.
capture_drawing <- function(expr) {
  drawn.lines <- list()
  drawn.polygons <- list()

  local_mocked_bindings(
    lines = function(x, y, ...) {
      drawn.lines[[length(drawn.lines) + 1]] <<- list(x = x, y = y)
      invisible(NULL)
    },
    polygon = function(x, y, ...) {
      drawn.polygons[[length(drawn.polygons) + 1]] <<- list(x = x, y = y)
      invisible(NULL)
    },
    .package = "EpiModel"
  )

  force(expr)
  list(lines = drawn.lines, polygons = drawn.polygons)
}

build_icm <- function() {
  set.seed(1)
  param <- param.icm(inf.prob = 0.2, act.rate = 0.5)
  init <- init.icm(s.num = 200, i.num = 5)
  control <- control.icm(type = "SI", nsteps = 10, nsims = 4, verbose = FALSE)
  icm(param, init, control)
}

build_netsim <- function() {
  set.seed(1)
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  control <- control.net(type = "SI", nsteps = 10, nsims = 4, verbose = FALSE)
  netsim(est, param.net(inf.prob = 0.3), init.net(i.num = 5), control)
}

test_that("sims restricts the mean line to the selected simulations", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (mod in list(build_icm(), build_netsim())) {
    all.sims <- rowMeans(mod$epi$i.num)
    two.sims <- rowMeans(mod$epi$i.num[, 1:2])
    expect_false(isTRUE(all.equal(all.sims, two.sims)))

    drawn <- capture_drawing(
      plot(mod, y = "i.num", sims = 1:2, mean.smooth = FALSE,
           sim.lines = FALSE, qnts = FALSE)
    )
    expect_length(drawn$lines, 1)
    expect_equal(drawn$lines[[1]]$y, two.sims)

    ## the default remains the mean over every simulation
    drawn <- capture_drawing(
      plot(mod, y = "i.num", mean.smooth = FALSE, sim.lines = FALSE,
           qnts = FALSE)
    )
    expect_equal(drawn$lines[[1]]$y, all.sims)

    ## a single simulation plots itself as the mean
    drawn <- capture_drawing(
      plot(mod, y = "i.num", sims = 3, mean.smooth = FALSE,
           sim.lines = FALSE, qnts = FALSE)
    )
    expect_equal(drawn$lines[[1]]$y, mod$epi$i.num[, 3])
  }
})

test_that("sims restricts the quantile band to the selected simulations", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (mod in list(build_icm(), build_netsim())) {
    expected <- apply(mod$epi$i.num[, 1:3], 1, quantile, probs = c(0.25, 0.75))

    drawn <- capture_drawing(
      plot(mod, y = "i.num", sims = 1:3, qnts = 0.5, qnts.smooth = FALSE,
           mean.line = FALSE, sim.lines = FALSE)
    )
    expect_length(drawn$polygons, 1)
    expect_equal(drawn$polygons[[1]]$y,
                 unname(c(expected[1, ], rev(expected[2, ]))))

    ## a single simulation has no band to draw
    drawn <- capture_drawing(
      plot(mod, y = "i.num", sims = 2, qnts = 0.5, mean.line = FALSE,
           sim.lines = FALSE)
    )
    expect_length(drawn$polygons, 0)
  }
})

test_that("sims restricts the simulation lines to the selected simulations", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (mod in list(build_icm(), build_netsim())) {
    drawn <- capture_drawing(
      plot(mod, y = "i.num", sims = c(2, 4), sim.lines = TRUE,
           mean.line = FALSE, qnts = FALSE)
    )
    expect_length(drawn$lines, 2)
    expect_equal(drawn$lines[[1]]$y, mod$epi$i.num[, 2])
    expect_equal(drawn$lines[[2]]$y, mod$epi$i.num[, 4])

    expect_error(plot(mod, y = "i.num", sims = 9), "Set sim to between 1 and 4")
    expect_error(plot(mod, y = "i.num", sims = 0), "Set sim to between 1 and 4")
  }
})

test_that("sims behaves the same through plot.epi.data.frame", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_netsim()
  df <- as.data.frame(mod)
  two.sims <- rowMeans(mod$epi$i.num[, 1:2])

  drawn <- capture_drawing(
    plot(df, y = "i.num", sims = 1:2, mean.smooth = FALSE, sim.lines = FALSE,
         qnts = FALSE)
  )
  expect_equal(drawn$lines[[1]]$y, two.sims)
})
