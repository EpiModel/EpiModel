context("multiple network models")

test_that("netsim runs with multiple networks", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

  dc10 <- dissolution_coefs(~offset(edges), 10, 0)
  dc20 <- dissolution_coefs(~offset(edges), 20, 0)
  dc1 <- dissolution_coefs(~offset(edges), 1, 0)
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  est10 <- est
  est1 <- update_dissolution(est, dc1)
  est20 <- update_dissolution(est, dc20)

  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  init <- init.net(i.num = 10)

  for (tergmLite in c(FALSE, TRUE)) {
    for (resimulate.network in unique(c(tergmLite, TRUE))) {
      control <- control.net(type = "SI", 
                             nsims = 3, 
                             nsteps = 5,
                             tergmLite = tergmLite,
                             resimulate.network = resimulate.network,
                             tergmLite.track.duration = TRUE,
                             nwstats.formula = multilayer(~triangle, "formation", ~mean.age, ~degree(0:3), "formation"),
                             verbose = TRUE)

      set.seed(0)
      mod <- netsim(list(est10, est1, est20, est1, est20), param, init, control)  

      stat_names <- list(c("triangle"),
                         c("edges", "nodematch.race"),
                         c("mean.age"),
                         c("degree0", "degree1", "degree2", "degree3"),
                         c("edges", "nodematch.race"))

      for (sim in seq_len(3)) {
        for (network in seq_len(5)) {
          expect_identical(colnames(get_nwstats(mod, sim = sim, network = network, mode = "list")[[1]]),
                           stat_names[[network]])
        }
      }

      expect_is(mod, "netsim")

      print(mod)

      plot(mod, type = "formation")

      if (tergmLite == FALSE) {
        plot(mod, type = "dissolution")
        plot(mod, type = "duration")
      }

      for (network in seq_len(5)) {
        print(mod, network = network)
        plot(mod, type = "formation", network = network)

        if (tergmLite == FALSE) {
          plot(mod, type = "dissolution", network = network)
          plot(mod, type = "duration", network = network)
        }
      }
    }
  }
})
