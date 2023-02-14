context("Extremal Targets")

test_that("extremal targets are handled correctly", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  nw <- set_vertex_attribute(nw, "sex", rbinom(50, 1, 0.5))

  target_stats <- c(25, 0, 10, 0, 0)
  target_stats_names <- c("edges", "nodematch.sex", "nodematch.race", "deg3+", "degree6")
  target_stats_names_regexp <- c("edges", "nodematch\\.sex", "nodematch\\.race", "deg3\\+", "degree6")

  est <- netest(nw, formation = ~edges + offset(degrange(3)) + nodematch("sex")
                                 + nodematch("race") + degrange(3) + offset(nodematch("sex"))
                                 + degree(6),
                target.stats = target_stats,
                coef.form = c(1, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_equal(est$target.stats, target_stats)
  expect_equal(est$target.stats.names, target_stats_names)

  dxs <- netdx(est, nsims = 10, dynamic = FALSE)
  print(dxs)
  plot(dxs)

  dxd <- netdx(est, nsims = 10, nsteps = 5, dynamic = TRUE)
  print(dxd)
  plot(dxd)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  print(mod)

  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "duration")
  plot(mod, type = "dissolution")
  plot(mod, type = "network")

  for (obj in list(dxs, dxd, mod)) {
    for(i in seq_along(target_stats_names_regexp)) {
      expect_output(print(obj), paste0(target_stats_names_regexp[i], " *", target_stats[i]))
    }
    expect_output(print(obj), paste0("offset\\(deg3\\+\\)", " *", "NA"))
    expect_output(print(obj), paste0("offset\\(nodematch\\.sex\\)", " *", "NA"))
  }
})
