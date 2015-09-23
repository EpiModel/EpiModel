context("Network diagnostics")

test_that("Edges only models", {
  num <- 50
  nw <- network.initialize(num, directed = FALSE)
  formation <- ~edges
  target.stats <- 15
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  ## Single simulation
  dx1 <- netdx(est1, nsims = 1, nsteps = 10, verbose = TRUE)
  expect_is(dx1, "netdx")
  print(dx1)
  plot(dx1)
  plot(dx1, method = "b")
  plot(dx1, type = "duration", mean.smooth = FALSE)
  plot(dx1, method = "b", type = "duration")
  plot(dx1, type = "dissolution")
  plot(dx1, method = "b", type = "dissolution")

  ## Multiple simulations
  dx2 <- netdx(est1, nsims = 2, nsteps = 10, verbose = TRUE)
  expect_is(dx2, "netdx")
  print(dx2)
  plot(dx2)
  plot(dx2, method = "b")
  plot(dx2, type = "duration")
  plot(dx2, method = "b", type = "duration")
  plot(dx2, type = "dissolution")
  plot(dx2, method = "b", type = "dissolution")

  ## Expanded monitoring formula
  dx3 <- netdx(est1, nsims = 2, nsteps = 10, verbose = TRUE,
               nwstats.formula = ~edges + concurrent)
  expect_is(dx3, "netdx")
  print(dx3)
  plot(dx3)
  plot(dx3, plots.joined = FALSE)
  plot(dx3, method = "b")
  plot(dx3, type = "duration")
  plot(dx3, method = "b", type = "duration")
  plot(dx3, type = "dissolution")
  plot(dx3, method = "b", type = "dissolution")

  ## Reduced monitoring formula
  dx4 <- netdx(est1, nsims = 2, nsteps = 10, verbose = TRUE,
               nwstats.formula = ~meandeg)
  expect_is(dx4, "netdx")
  print(dx4)
  plot(dx4)
  plot(dx4, method = "b")
  plot(dx4, type = "duration")
  plot(dx4, method = "b", type = "duration")
  plot(dx4, type = "dissolution")
  plot(dx4, method = "b", type = "dissolution")
})


test_that("Offset terms", {
  n <- 50
  nw <- network.initialize(n, directed = FALSE)
  nw <- set.vertex.attribute(nw, "loc", rep(0:1, each = n/2))
  dissolution <- ~offset(edges)
  duration <- 40
  coef.diss <- dissolution_coefs(dissolution, duration)
  formation <- ~edges + offset(nodemix("loc", base = c(1, 3)))
  target.stats <- 15
  est2 <- netest(nw, formation, target.stats, coef.diss, coef.form = -Inf,
                 verbose = FALSE)
  dx <- netdx(est2, nsims = 2, nsteps = 50, verbose = TRUE)
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
  dx <- netdx(est2, nsims = 2, nsteps = 50, verbose = TRUE,
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
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 40)
  formation <- ~edges + nodemix("loc", base = c(1, 3))
  target.stats <- c(15, 0)
  est3 <- netest(nw, formation, target.stats, coef.diss, verbose = TRUE)

  dx <- netdx(est3, nsims = 2, nsteps = 50, verbose = FALSE)
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


test_that("More complicated faux offset term", {
  # skip_on_cran()
  nw <- network.initialize(1000, directed = FALSE)
  nw <- set.vertex.attribute(nw, "sexor",
                             sample(c(rep(1,20), rep(2,460),
                                      rep(3,20), rep(4,500))))
  nw <- set.vertex.attribute(nw, "region", sample(rep(1:5,200)))
  fit <- netest(nw,
                formation = ~edges +
                             nodemix("sexor", base = 1) +
                             nodematch("region"),
                target.stats = c(463, 0, 0, 18, 0, 6, 0, 380, 25, 0, 400),
                coef.diss = dissolution_coefs(~offset(edges), 60))
  dx <- netdx(fit, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
})


test_that("Static diagnostic simulations", {
  nw <- network.initialize(100, directed = FALSE)
  formation <- ~edges + concurrent
  target.stats <- c(50, 20)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est4 <- netest(nw, formation, target.stats, coef.diss, verbose = TRUE)

  dx <- netdx(est4, dynamic = FALSE, nsims = 250,
              nwstats.formula = ~edges + meandeg + concurrent)
  expect_is(dx, "netdx")

  plot(dx)
  plot(dx, stats = "meandeg")
  plot(dx, plots.joined = FALSE)

  plot(dx, method = "b", col = "bisque")

  expect_error(plot(dx, method = "b", type = "duration"))
  expect_error(plot(dx, method = "b", type = "dissolution"))
})


test_that("Parallel methods", {
  skip_on_cran()
  num <- 50
  nw <- network.initialize(num, directed = FALSE)
  formation <- ~edges
  target.stats <- 15
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  ## Single simulation
  dx1 <- netdx(est1, nsims = 2, nsteps = 25, ncores = 2)
  expect_is(dx1, "netdx")
})
