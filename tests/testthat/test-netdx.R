context("Network diagnostics (All SOC)")

for (trim in c(FALSE, TRUE)) {

  test_that("Edges only models", {
    skip_on_cran()
    num <- 50
    nw <- network_initialize(n = num)
    formation <- ~edges
    target.stats <- 15
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
      est1 <- trim_netest(est1)
    }
    ## Single simulation
    dx1 <- netdx(est1, nsims = 1, nsteps = 10, verbose = FALSE)
    expect_is(dx1, "netdx")
    capture_output(
      print(dx1)
    )
    plot(dx1)
    plot(dx1, method = "b")
    plot(dx1, type = "duration", mean.smooth = FALSE)
    plot(dx1, method = "b", type = "duration")
    plot(dx1, type = "dissolution")
    plot(dx1, method = "b", type = "dissolution")

    ## Multiple simulations
    dx2 <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE)
    expect_is(dx2, "netdx")
    capture_output(
      print(dx2)
    )
    plot(dx2)
    plot(dx2, method = "b")
    plot(dx2, type = "duration")
    plot(dx2, method = "b", type = "duration")
    plot(dx2, type = "dissolution")
    plot(dx2, method = "b", type = "dissolution")

    ## Expanded monitoring formula
    dx3 <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE,
                 nwstats.formula = ~edges + concurrent)
    expect_is(dx3, "netdx")
    capture_output(
      print(dx3)
    )
    plot(dx3)
    plot(dx3, plots.joined = FALSE)
    plot(dx3, method = "b")
    plot(dx3, type = "duration")
    plot(dx3, method = "b", type = "duration")
    plot(dx3, type = "dissolution")
    plot(dx3, method = "b", type = "dissolution")
    plot(dx3, type = "formation", sim.lines = FALSE, plots.joined = FALSE,
         mean.line = FALSE, qnts = FALSE)

    ## Reduced monitoring formula
    dx4 <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE,
                 nwstats.formula = ~meandeg)
    expect_is(dx4, "netdx")
    capture_output(
      print(dx4)
    )
    plot(dx4)
    plot(dx4, method = "b")
    plot(dx4, type = "duration")
    plot(dx4, method = "b", type = "duration")
    plot(dx4, type = "dissolution")
    plot(dx4, method = "b", type = "dissolution")
  })

  test_that("Formation plot color vector length", {
    skip_on_cran()
    n <-  100
    mean.degree <- ((0 * 0.10) + (1 * 0.41) + (2 * 0.25) + (3 * 0.22))
    expected.concurrent <- n * 0.49
    expected.edges <- (mean.degree) * (n / 2)
    nw <- network_initialize(n = n)
    formation <- ~edges + concurrent + degrange(from = 4)
    target.stats <- c(expected.edges, expected.concurrent, 0)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 40,
                                   d.rate = 0.01)
    est <- netest(nw, formation, target.stats, coef.diss)
    if (trim == TRUE) {
     est <- trim_netest(est)
    }
    dx <- netdx(est, nsims = 2, nsteps = 500,
                nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent, verbose = FALSE)

    expect_error(plot(dx, sim.col = c("green", "orange")))
    expect_error(plot(dx, mean.col = c("green", "orange")))
    expect_error(plot(dx, targ.col = c("green", "orange")))
    expect_error(plot(dx, qnts.col = c("green", "orange")))
  })

  test_that("Netdx duration and dissolution plots error when
            skip.dissolution = TRUE", {
    skip_on_cran()
    nw <- network_initialize(n = 100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2)
    est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est <- trim_netest(est)
    }
    dx2 <- netdx(est, nsims = 1, nsteps = 500, skip.dissolution = TRUE, verbose = FALSE)

    expect_error(plot(dx2, type = "duration"),
                 "Plots of type duration and dissolution only available if netdx run with skip.dissolution = FALSE")
    expect_error(plot(dx2, type = "dissolution"),
                 "Plots of type duration and dissolution only available if netdx run with skip.dissolution = FALSE")
  })

  test_that("Offset terms", {
    skip_on_cran()
    n <- 50
    nw <- network_initialize(n = n)
    nw <- set_vertex_attribute(nw, "loc", rep(0:1, each = n / 2))
    dissolution <- ~offset(edges)
    duration <- 40
    coef.diss <- dissolution_coefs(dissolution, duration)
    formation <- ~edges + offset(nodemix("loc", base = c(1, 3)))
    target.stats <- 15
    est2 <- netest(nw, formation, target.stats, coef.diss, coef.form = -Inf,
                   verbose = FALSE)
    if (trim == TRUE) {
     est2 <- trim_netest(est2)
    }
    dx <- netdx(est2, nsims = 2, nsteps = 50, verbose = FALSE)
    expect_is(dx, "netdx")
    capture_output(
      print(dx)
    )
    plot(dx)
    plot(dx, plots.joined = FALSE)
    plot(dx, method = "b")
    plot(dx, type = "duration")
    plot(dx, method = "b", type = "duration")
    plot(dx, type = "dissolution")
    plot(dx, method = "b", type = "dissolution")
    rm(dx)

    ## Offset term with expanded formula
    dx <- netdx(est2, nsims = 2, nsteps = 50, verbose = FALSE,
                nwstats.formula = ~edges + meandeg + concurrent + nodematch("loc"))
    expect_is(dx, "netdx")
    capture_output(
      print(dx)
    )
    plot(dx)
    plot(dx, plots.joined = FALSE)
    plot(dx, method = "b")
    plot(dx, type = "duration")
    plot(dx, method = "b", type = "duration")
    plot(dx, type = "dissolution")
    plot(dx, method = "b", type = "dissolution")
  })

  test_that("Faux offset term", {
    skip_on_cran()
    n <- 50
    nw <- network_initialize(n = n)
    nw <- set_vertex_attribute(nw, "loc", rep(0:1, each = n / 2))
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 40)
    formation <- ~edges + nodemix("loc", base = c(1, 3))
    target.stats <- c(15, 0)
    est3 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est3 <- trim_netest(est3)
    }

    dx <- netdx(est3, nsims = 2, nsteps = 50, verbose = FALSE)
    expect_is(dx, "netdx")
    capture_output(
      print(dx)
    )
    plot(dx)
    plot(dx, plots.joined = FALSE)
    plot(dx, method = "b")
    plot(dx, type = "duration")
    plot(dx, method = "b", type = "duration")
    plot(dx, type = "dissolution")
    plot(dx, method = "b", type = "dissolution")
  })


  test_that("More complicated faux offset term", {
    skip_on_cran()
    nw <- network_initialize(n = 1000)
    nw <- set_vertex_attribute(nw, "sexor",
                               sample(c(rep(1, 20), rep(2, 460),
                                        rep(3, 20), rep(4, 500))))
    nw <- set_vertex_attribute(nw, "region", sample(rep(1:5, 200)))
    fit <- netest(nw,
                  formation =
                    ~edges + nodemix("sexor", levels2 = c(-1, -3, -5)) + nodematch("region"),
                  target.stats = c(463, 0, 0, 0, 10, 210, 10, 110, 90),
                  coef.diss = dissolution_coefs(~offset(edges), 60),
                  set.control.ergm = control.ergm(MCMLE.termination = "Hummel",
                                                  MCMLE.samplesize = 1024,
                                                  MCMLE.interval = 1024,
                                                  MCMLE.effectiveSize = NULL))
    if (trim == TRUE) {
     fit <- trim_netest(fit)
    }
    dx <- netdx(fit, nsteps = 10, verbose = FALSE)
    capture_output(
      print(dx)
    )
    expect_is(dx, "netdx")
  })


  test_that("Static diagnostic simulations", {
    skip_on_cran()
    nw <- network_initialize(n = 100)
    formation <- ~edges + concurrent
    target.stats <- c(50, 20)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
    est4 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est4 <- trim_netest(est4)
    }

    dx <- netdx(est4, dynamic = FALSE, nsims = 250,
                nwstats.formula = ~edges + meandeg + concurrent, verbose = FALSE)
    expect_is(dx, "netdx")

    plot(dx)
    plot(dx, stats = "meandeg")
    plot(dx, plots.joined = FALSE)

    plot(dx, method = "b", col = "bisque")

    expect_error(plot(dx, method = "b", type = "duration"))
    expect_error(plot(dx, method = "b", type = "dissolution"))

    # test for default formation model
    dx <- netdx(est4, dynamic = FALSE, nsims = 250, verbose = FALSE)
    expect_is(dx, "netdx")
  })


  test_that("Parallel methods", {
    skip_on_cran()
    skip_on_os(os = "windows")
    num <- 50
    nw <- network_initialize(n = num)
    formation <- ~edges
    target.stats <- 15
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est1 <- trim_netest(est1)
    }

    dx1 <- netdx(est1, nsims = 2, nsteps = 25, ncores = 2, verbose = FALSE)
    expect_is(dx1, "netdx")
  })

  test_that("error checking", {
    skip_on_cran()
    nw <- network_initialize(n = 25)
    est <- netest(nw, formation = ~edges, target.stats = 25,
                  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                  edapprox = TRUE, verbose = FALSE)
    if (trim == TRUE) {
     est <- trim_netest(est)
    }
    expect_error(netdx(x = 1, nsteps = 100))
    expect_error(netdx(est), "Specify number of time steps with nsteps")
  })

  test_that("Cross sectional ergm dynamic error check", {
    skip_on_cran()
    nw <- network_initialize(n = 100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
    est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est <- trim_netest(est)
    }
    expect_error(netdx(est, nsims = 5, nsteps = 500),
                 "Running dynamic diagnostics on a cross-sectional")
  })

  test_that("print.netdx output", {
    skip_on_cran()
    nw <- network_initialize(n = 100)
    formation <- ~edges + concurrent
    target.stats <- c(50, 25)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
    est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
     est <- trim_netest(est)
    }
    dx <- netdx(est, nsims = 1, nsteps = 10, verbose = FALSE)
    expect_output(print(dx), "edges")
    expect_output(print(dx), "concurrent")
    expect_output(print(dx), "Duration Diagnostics")
    expect_output(print(dx), "Dissolution Diagnostics")
    expect_output(print(dx), "Target")
    expect_output(print(dx), "Sim Mean")
    expect_output(print(dx), "Pct Diff")
    expect_output(print(dx), "Sim SE")

    dx <- netdx(est, nsims = 1, nsteps = 100, verbose = FALSE)
    expect_output(print(dx), "NA")
  })

  test_that("print.netdx and plot.netdx with heterogeneous diss", {
    skip_on_cran()
    set.seed(12345)
    nw <- network_initialize(n = 100)
    nw <- set_vertex_attribute(nw, "neighborhood", rep(1:10, 10))
    nw <- set_vertex_attribute(nw, "position", rep(1:4, 25))
    formation <- ~edges+nodematch("neighborhood", diff = TRUE)
    target.stats <- c(100, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
    coef.diss <- dissolution_coefs(dissolution =
                   ~offset(edges)+offset(nodematch("neighborhood", diff = TRUE)),
                     duration = c(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30))
    est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
      est <- trim_netest(est)
    }
    dx11 <- netdx(est, nsims = 5, nsteps = 100, verbose = FALSE)
    expect_length(dx11$stats.table.duration$Target, 11)
    expect_length(dx11$stats.table.dissolution$`Sim Mean`, 11)
    suppressWarnings({
      expect_output(print(dx11), "match.neighborhood.7")
      plot(dx11)
      plot(dx11, type = "duration")
      plot(dx11, type = "dissolution")
    })

    formation <- ~edges + nodemix("position")
    target.stats <- c(100, 9, 7, 8, 8, 12, 8, 15, 5, 10)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodemix("position")),
                                   duration = c(80, 75, 80, 70, 85, 65, 70, 75, 60, 80))
    est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    if (trim == TRUE) {
      est <- trim_netest(est)
    }
    dx12 <- netdx(est, nsims = 7, nsteps = 30, verbose = FALSE)
    expect_output(print(dx12), "60")
    expect_output(print(dx12, digits = 4), "0.0133")
    expect_output(print(dx12), "mix.position.1.4")
    plot(dx12, qnts = FALSE)
    plot(dx12, type = "duration", qnts=0.8)
    plot(dx12, type = "dissolution", sim.col = 1:10, verbose = FALSE)

    dx13 <- netdx(est, nsims = 1, nsteps = 1, verbose = FALSE)
    expect_true(dim(dx13$pages)[1]==1 & dim(dx13$pages)[2]==10 & dim(dx13$pages)[3]==1)
    expect_true(dim(dx13$prop.diss)[1]==1 & dim(dx13$prop.diss)[2]==10 & dim(dx13$prop.diss)[3]==1)
  })
}

test_that("z scores are not large for a reasonably long simulation", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rep(0:1, length.out = 100))
  est <- netest(nw, formation = ~edges + nodematch("race", diff = TRUE),
                target.stats = c(50, 10, 10),
                coef.diss = dissolution_coefs(~offset(edges) +
                                                offset(nodematch("race", diff = TRUE)),
                                              c(10, 20, 15)),
                verbose = FALSE
  )

  dx <- netdx(est, nsteps = 1000, nsims = 1, verbose = FALSE)

  expect_true(all(abs(dx$stats.table.formation[["Z Score"]]) < 20))
  expect_true(all(abs(dx$stats.table.duration[["Z Score"]]) < 20))
  expect_true(all(abs(dx$stats.table.dissolution[["Z Score"]]) < 20))

  dxs <- netdx(est, nsims = 1000, dynamic = FALSE, verbose = FALSE)

  expect_true(all(abs(dxs$stats.table.formation[["Z Score"]]) < 20))
})

test_that("make_stats_table behaves as expected", {
  skip_on_cran()
  stat_names <- letters[1:5]
  stats_1 <- 1:5
  stats_2 <- c(6, 9, 5, 2, 7)
  m1 <- matrix(stats_1, ncol = 5, nrow = 10, byrow = TRUE)
  m2 <- matrix(stats_2, ncol = 5, nrow = 10, byrow = TRUE)
  colnames(m2) <- colnames(m1) <- stat_names
  stats <- list(m1, m2)
  targets <- c(a = 3, f = 6, e = 7, c = 9, x = 10, d = 2000)
  targs <- c(3, NA, 9, 2000, 7)
  mst <- make_stats_table(stats, targets)
  expect_equal(rownames(mst), colnames(m1))
  expect_equal(mst[["Sim Mean"]], (stats_1 + stats_2)/2)
  expect_equal(mst[["Target"]], targs)
  expect_equal(mst[["Pct Diff"]], 100*((stats_1 + stats_2)/2 - targs)/targs)
})
