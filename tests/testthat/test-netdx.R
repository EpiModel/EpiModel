context("Network diagnostics")

test_that("Edges only models", {

  num <- 50
  nw <- network.initialize(num, directed = FALSE)
  formation <- ~ edges
  target.stats <- 15
  dissolution <- ~ offset(edges)
  duration <- 20
  coef.diss <- dissolution_coefs(dissolution, duration)
  est1 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 verbose = FALSE)

  ## Single simulation
  dx <- netdx(est1, nsims = 1, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")
  rm(dx)

  ## Multiple simulations
  dx <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")
  rm(dx)

  ## Expanded monitoring formula
  dx <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE,
              nwstats.formula = ~edges + concurrent)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, plots.joined = FALSE)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")
  rm(dx)

  ## Reduced monitoring formula
  dx <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE,
              nwstats.formula = ~meandeg)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")

})


test_that("Offset terms", {

  n <- 50
  nw <- network.initialize(n, directed = FALSE)
  nw <- set.vertex.attribute(nw, "loc", rep(0:1, each = n/2))
  dissolution <- ~ offset(edges)
  duration <- 40
  coef.diss <- dissolution_coefs(dissolution, duration)
  formation <- ~ edges + offset(nodemix("loc", base=c(1, 3)))
  target.stats <- 15
  est2 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 coef.form = -Inf,
                 verbose = FALSE)

  dx <- netdx(est2, nsims = 2, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, plots.joined = FALSE)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")
  rm(dx)

  ## Offset term with expanded formula
  dx <- netdx(est2, nsims = 2, nsteps = 10, verbose = FALSE,
              nwstats.formula = ~edges + meandeg + concurrent + nodematch("loc"))
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, plots.joined = FALSE)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")

})


test_that("Faux offset term", {

  n <- 50
  nw <- network.initialize(n, directed = FALSE)
  nw <- set.vertex.attribute(nw, "loc", rep(0:1, each = n/2))
  dissolution <- ~ offset(edges)
  duration <- 40
  coef.diss <- dissolution_coefs(dissolution, duration)
  formation <- ~ edges + nodemix("loc", base=c(1, 3))
  target.stats <- c(15, 0)
  est3 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 verbose = FALSE)

  dx <- netdx(est3, nsims = 2, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
  print(dx)
  plot(dx)
  plot(dx, plots.joined = FALSE)
  plot(dx, method = "b")
  plot(dx, type = "duration")
  plot(dx, method = "b", type = "duration")
  plot(dx, type = "dissolution")
  plot(dx, method = "b", type = "dissolution")

})


test_that("Static diagnostic simulations", {

  nw <- network.initialize(100, directed = FALSE)
  formation <- ~ edges + concurrent
  dissolution <- ~ offset(edges)
  target.stats <- c(50, 20)
  coef.diss <- dissolution_coefs(dissolution, duration = 10)

  est4 <- netest(nw, formation, dissolution,
                 target.stats, coef.diss, verbose = FALSE)

  dx <- netdx(est4, dynamic = FALSE, nsims = 250,
              nwstats.formula = ~ edges + meandeg + concurrent)
  expect_is(dx, "netdx")

  plot(dx)
  plot(dx, stats = "meandeg")
  plot(dx, plots.joined = FALSE)

  plot(dx, method = "b", col = "bisque")

  expect_error(plot(dx, method = "b", type = "duration"))
  expect_error(plot(dx, method = "b", type = "dissolution"))

})


