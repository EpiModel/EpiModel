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

  est <- list(est10, est1, est20, est1, est20)

  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  init <- init.net(i.num = 10)

  nsims <- 3L
  nnets <- 5L
  tergm_nets <- c(1,3,5)

  stat_names <- list(c("triangle"),
                     c("edges", "nodematch.race"),
                     c("mean.age"),
                     c("degree0", "degree1", "degree2", "degree3"),
                     c("edges", "nodematch.race"))

  for (tergmLite in c(FALSE, TRUE)) {
    for (resimulate.network in unique(c(tergmLite, TRUE))) {
      if (tergmLite == TRUE) {
        save.other <- c("nw", "attr", "el", "temp")
      } else {
        save.other <- c("nw", "attr", "temp")      
      }

      for (iteration in 1:2) {
        if (iteration == 1) {
          nsteps <- 5L
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
                                 save.other = save.other)
          print(control)
          basis <- est
        } else {
          if (tergmLite == TRUE) {
            next
          }
          control$start <- nsteps + 1L
          nsteps <- 11L
          control$nsteps <- nsteps
          basis <- sim          
        }
        sim <- netsim(basis, param, init, control)

        print(sim)
        plot(sim)
        plot(sim, type = "formation")
        if (tergmLite == FALSE) {
          plot(sim, type = "dissolution")
          plot(sim, type = "duration")
        }

        for (network in seq_len(nnets)) {
          print(sim, network = network)
          plot(sim, network = network)
          plot(sim, type = "formation", network = network)
          if (tergmLite == FALSE) {
            plot(sim, type = "dissolution", network = network)
            plot(sim, type = "duration", network = network)
          }
        }

        for(simno in seq_len(nsims)) {
          for(network in seq_len(nnets)) {
            if (tergmLite == TRUE) {
              expect_is(sim$nw[[simno]][[network]], "networkLite")
              expect_equal(sim$el[[simno]][[network]], as.edgelist(sim$nw[[simno]][[network]]))
              if (network %in% tergm_nets) {
                expect_equal(sim$nw[[simno]][[network]] %n% "time", nsteps)
                expect_true(all((sim$nw[[simno]][[network]] %n% "lasttoggle")[,3] >= 0))
                expect_true(all((sim$nw[[simno]][[network]] %n% "lasttoggle")[,3] <= nsteps))
              }
            } else {
              expect_is(sim$nw[[simno]][[network]], "networkDynamic")            
            }

            expect_equal(sim$nwparam[[network]]$coef.form[1],
                         est[[network]]$coef.form[1] +
                           log(network.size(est[[network]]$newnetwork)) -
                           log(network.size(sim$nw[[simno]][[network]])))

            expect_equal(network.size(est[[network]]$newnetwork), sim$epi$num[[simno]][1])
            expect_equal(network.size(sim$nw[[simno]][[network]]), sim$epi$num[[simno]][nsteps])

            stats_matrix <- get_nwstats(sim, network = network, mode = "list")[[simno]]
            expect_equal(NCOL(stats_matrix), length(stat_names[[network]]))
            expect_equal(NROW(stats_matrix), nsteps)
            expect_equal(colnames(stats_matrix), stat_names[[network]])

            final_stats <- stats_matrix[nsteps,]
            names(final_stats) <- colnames(stats_matrix)
            if (tergmLite == TRUE) {
              nwL <- networkLite(sim$el[[simno]][[network]], sim$attr[[simno]])
              if (network %in% tergm_nets) {
                nwL %n% "time" <- sim$nw[[simno]][[network]] %n% "time"
                nwL %n% "lasttoggle" <- sim$nw[[simno]][[network]] %n% "lasttoggle"
              }
              expect_equal(final_stats,
                           summary(get_network_control(sim, network = network, "nwstats.formula"), basis = nwL))
            } else {
              expect_equal(final_stats,
                           summary(get_network_control(sim, network = network, "nwstats.formula"), at = nsteps, basis = sim$nw[[simno]][[network]])[,,drop=TRUE])
            }
          }
        }
      }
    }
  }
})
