context("multiple network models")

test_that("netsim runs with multiple networks, with open or closed population", {
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

  param_open <- param.net(inf.prob = 0.3, act.rate = 0.1, a.rate = 0.03,
                          di.rate = 0.03, ds.rate = 0.03)
  param_closed <- param.net(inf.prob = 0.3, act.rate = 0.1)

  init <- init.net(i.num = 10)

  nsims <- 3L
  nnets <- 5L
  tergm_nets <- c(1,3,5)

  resim_stat_names <- list(c("triangle", "nodefactor.deg.4.1"),
                           c("edges", "nodematch.race"),
                           c("mean.age", "nodefactor.deg.2.2"),
                           c("degree0", "degree1", "degree2", "degree3"),
                           c("edges", "nodematch.race", "nodefactor.deg.5.1"))

  resim_nwstats_formulas <- multilayer(~triangle + nodefactor("deg.4", levels = I(1)),
                                       "formation",
                                       ~mean.age + nodefactor("deg.2", levels = I(2)),
                                       ~degree(0:3),
                                       ~edges + nodematch("race") + nodefactor("deg.5", levels = I(1)))

  stat_names <- list(c("triangle"),
                     c("edges", "nodematch.race"),
                     c("mean.age"),
                     c("degree0", "degree1", "degree2", "degree3"),
                     c("edges", "nodematch.race"))

  nwstats_formulas <- multilayer(~triangle,
                                 "formation",
                                 ~mean.age,
                                 ~degree(0:3),
                                 "formation")

  for (tergmLite in c(FALSE, TRUE)) {
    for (resimulate.network in unique(c(tergmLite, TRUE))) {
      for (open_population in unique(c(FALSE, resimulate.network))) {
        if (tergmLite == TRUE) {
          save.other <- c()
        } else {
          save.other <- c()
        }
        if (open_population == TRUE) {
          param <- param_open
        } else {
          param <- param_closed
        }
        if (resimulate.network == TRUE) {
          sim_stat_names <- resim_stat_names
          sim_nwstats_formulas <- resim_nwstats_formulas
        } else {
          sim_stat_names <- stat_names
          sim_nwstats_formulas <- nwstats_formulas
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
                                   dat.updates = function (dat, at, network) {
                                     if (network > 0L) {
                                       if (get_control(dat, "tergmLite") == TRUE) {
                                         dat <- set_attr(dat, paste0("deg.", network),
                                                         get_degree(dat$run$el[[network]]))
                                       } else {
                                         deg_attr <- get_degree(as.edgelist(network.collapse(dat$nw[[network]],
                                                                                             at = at,
                                                                                             retain.all.vertices = TRUE)))
                                         for (other_net in seq_len(dat$num.nw)) {
                                           dat$nw[[other_net]] <- set_vertex_attribute(dat$nw[[other_net]],
                                                                                       paste0("deg.", network),
                                                                                       deg_attr)
                                         }
                                       }
                                     }
                                     dat
                                   },
                                   nwstats.formula = sim_nwstats_formulas,
                                   verbose = TRUE,
                                   save.network = TRUE,
                                   save.run = TRUE,
                                   save.other = c(save.other))
            print(control)
            basis <- est
          } else {
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
            expect_equal(sim$nwparam[[network]]$coef.form[1],
                         est[[network]]$coef.form[1] +
                           log(network.size(est[[network]]$newnetwork)) -
                           log(sim$run[[1]]$num),
                         tolerance = 1e-6)
          }

          for(simno in seq_len(nsims)) {
            for(network in seq_len(nnets)) {
              expect_equal(sim$coef.form[[simno]][[network]][1],
                           est[[network]]$coef.form[1] +
                             log(network.size(est[[network]]$newnetwork)) -
                             log(sim$run[[simno]]$num),
                           tolerance = 1e-6)
              if (tergmLite == TRUE) {
                expect_is(sim$network[[simno]][[network]], "networkLite")
                if (open_population == FALSE) {
                  expect_equal(sim$run[[simno]]$el[[network]], as.edgelist(sim$network[[simno]][[network]]))
                  expect_equal(sim$run[[simno]]$attr[[paste0("deg.", network)]], get_degree(sim$run[[simno]]$el[[network]]))
                }
                if (network %in% tergm_nets) {
                  expect_equal(sim$network[[simno]][[network]] %n% "time", nsteps)
                  expect_true(all((sim$network[[simno]][[network]] %n% "lasttoggle")[,3] >= 0))
                  expect_true(all((sim$network[[simno]][[network]] %n% "lasttoggle")[,3] <= nsteps))
                }
              } else {
                expect_is(sim$network[[simno]][[network]], "networkDynamic")
                if (open_population == FALSE) {
                  if (resimulate.network == TRUE || iteration == 2) {
                    expect_equal(sim$network[[simno]][[network]] %v% paste0("deg.", network),
                                 get_degree(as.edgelist(network.collapse(sim$network[[simno]][[network]], at = nsteps))))
                  }
                }
              }

              expect_equal(network.size(est[[network]]$newnetwork), sim$epi$num[[simno]][1])
              if (tergmLite == TRUE) {
                if (open_population == FALSE) {
                  expect_equal(sim$epi$num[nsteps - 1,simno],
                               sim$run[[simno]]$num)
                }
              } else {
                expect_equal(network.size(network.collapse(sim$network[[simno]][[network]], at = nsteps - 1)),
                               sim$run[[simno]]$num)
                expect_equal(network.size(network.collapse(sim$network[[simno]][[network]], at = nsteps - 1)),
                               sim$run[[simno]]$num)
              }

              stats_matrix <- get_nwstats(sim, network = network, mode = "list")[[simno]]
              expect_equal(NCOL(stats_matrix), length(sim_stat_names[[network]]))
              expect_equal(NROW(stats_matrix), nsteps)
              expect_equal(colnames(stats_matrix), sim_stat_names[[network]])

              if (open_population == FALSE) {
                final_stats <- stats_matrix[nsteps,]
                names(final_stats) <- colnames(stats_matrix)
                if (tergmLite == TRUE) {
                  nwL <- networkLite(sim$run[[simno]]$el[[network]], sim$run[[simno]]$attr)
                  if (network %in% tergm_nets) {
                    nwL %n% "time" <- sim$network[[simno]][[network]] %n% "time"
                    nwL %n% "lasttoggle" <- sim$network[[simno]][[network]] %n% "lasttoggle"
                  }
                  expect_equal(final_stats,
                               summary(get_network_control(sim, network = network, "nwstats.formula"), basis = nwL))
                } else {
                  expect_equal(final_stats,
                               summary(get_network_control(sim, network = network, "nwstats.formula"), at = nsteps, basis = sim$network[[simno]][[network]])[,,drop=TRUE])
                }
              }
            }
          }
        }
      }
    }
  }
})

test_that("multilayer specifications", {
  nw <- network_initialize(n = 100)

  dc1 <- dissolution_coefs(~offset(edges), 1, 0)
  dc100 <- dissolution_coefs(~offset(edges), 100, 0)
  dc200 <- dissolution_coefs(~offset(edges), 200, 0)

  est1 <- netest(nw, formation = ~edges,
                 target.stats = c(250),
                 coef.diss = dc1,
                 verbose = FALSE)

  est100 <- update_dissolution(est1, dc100)
  est200 <- update_dissolution(est1, dc200)

  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI",
                         nsims = 2,
                         nsteps = 10,
                         tergmLite = TRUE,
                         resimulate.network = TRUE,
                         tergmLite.track.duration = multilayer(TRUE, FALSE, TRUE, FALSE, FALSE),
                         nwstats.formula = multilayer(~edges + triangle + mean.age,
                                                      ~edges + degree(0:1) + mean.age,
                                                      ~edges + degree(1:2) + mean.age,
                                                      ~edges + concurrent + mean.age,
                                                      "formation"),
                         verbose = FALSE)

  basis <- list(est100, est200, est1, est1, est200)
  sim <- netsim(basis, param, init, control)

  for (network in seq_len(5L)) {
    stats <- get_nwstats(sim, network = network)
    if (network < 5L) {
      expect_true(all(stats$mean.age >= 0))
      if (get_network_control(sim, network, "tergmLite.track.duration") == TRUE) {
        expect_true(all(stats$mean.age <= 11))
      } else {
        expect_true(all(stats$mean.age >= 1000)) # large number indicating tergm defaults are being used
      }
    }
  }
})
