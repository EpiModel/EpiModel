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

  nsims <- 3L
  nsteps <- 5L

  for (tergmLite in c(FALSE, TRUE)) {
    for (resimulate.network in unique(c(tergmLite, TRUE))) {
      control <- control.net(type = "SI", 
                             nsims = nsims, 
                             nsteps = nsteps,
                             tergmLite = tergmLite,
                             resimulate.network = resimulate.network,
                             tergmLite.track.duration = TRUE,
                             dat.updates = function (dat, at, network) dat,
                             nwstats.formula = multilayer(~triangle, "formation", ~mean.age, ~degree(0:3), "formation"),
                             verbose = TRUE,
                             save.network = !tergmLite,
                             save.other = c("attr"))
      print(control)

      mod <- netsim(list(est10, est1, est20, est1, est20), param, init, control)  

      stat_names <- list(c("triangle"),
                         c("edges", "nodematch.race"),
                         c("mean.age"),
                         c("degree0", "degree1", "degree2", "degree3"),
                         c("edges", "nodematch.race"))

      for (sim in seq_len(3)) {
        for (network in seq_len(5)) {
          stats_matrix <- get_nwstats(mod, sim = sim, network = network, mode = "list")[[1]]
          expect_identical(NCOL(stats_matrix), length(stat_names[[network]]))
          expect_identical(NROW(stats_matrix), nsteps)
          expect_identical(colnames(stats_matrix), stat_names[[network]])
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
      
      if (tergmLite == FALSE) {
        control$start <- 6L
        control$nsteps <- 11L

        print(control)

        mod2 <- netsim(mod, param, init, control)  

        for (sim in seq_len(3)) {
          for (network in seq_len(5)) {
            stats_matrix <- get_nwstats(mod2, sim = sim, network = network, mode = "list")[[1]]
            expect_identical(NCOL(stats_matrix), length(stat_names[[network]]))
            expect_identical(NROW(stats_matrix), control$nsteps)
            expect_identical(colnames(stats_matrix), stat_names[[network]])
          }
        }

        expect_is(mod2, "netsim")
        
        print(mod2)
        
        plot(mod2, type = "formation")        
        plot(mod2, type = "dissolution")
        plot(mod2, type = "duration")
        
        for (network in seq_len(5)) {
          print(mod, network = network)
          plot(mod, type = "formation", network = network)        
          plot(mod, type = "dissolution", network = network)
          plot(mod, type = "duration", network = network)
        }
      }
    }
  }
})
