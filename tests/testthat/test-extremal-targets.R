context("Extremal Targets")

test_that("extremal targets are handled correctly", {
  ff <- ~edges + offset(degrange(6)) + nodematch("sex") + nodematch("race") +
          degrange(6) + offset(nodematch("sex")) + degree(5)

  ff_sum <- ~edges + nodematch("sex") + nodematch("race") + degrange(6) +
              degree(5)

  for (egor in list(FALSE, TRUE)) {
    nw <- network_initialize(n = 100)
    nw <- set_vertex_attribute(nw, "race", rbinom(100, 1, 0.5))
    nw <- set_vertex_attribute(nw, "sex", rbinom(100, 1, 0.5))

    target_stats <- c(50, 0, 20, 0, 0)
    if (egor == TRUE) {
      nw <- san(ff, basis = nw, target.stats = target_stats, offset.coef = c(-Inf,-Inf))
      target_stats <- summary(ff_sum, basis = nw)
      nw <- ergm.ego:::as.egor.network(nw)
    }

    target_stats_names <- c("edges", "nodematch.sex", "nodematch.race", "deg6+", "degree5")
    target_stats_names_regexp <- c("edges", "nodematch\\.sex", "nodematch\\.race", "deg6\\+", "degree5")

    est <- netest(nw, formation = ff,
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
      expect_output(print(obj), paste0("offset\\(deg6\\+\\)", " *", "NA"))
      expect_output(print(obj), paste0("offset\\(nodematch\\.sex\\)", " *", "NA"))
    }
  }
})
