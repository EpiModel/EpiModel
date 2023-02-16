context("ergm.ego (All SOC)")

# Test ergm.ego ---------------------------------------------------

test_that("ergm and ergm.ego produce the same results in EpiModel", {
  skip_on_cran()
  for (trim in list(FALSE, TRUE)) {
    for (ppopnw in list(FALSE, TRUE)) {
      nw <- network_initialize(n = 50)
      nw %v% "race" <- rep(1:2, length.out = 50)
      nw <- san(nw ~ edges, target.stats = c(50))

      nwe <- as.egor(nw)

      nw[,] <- FALSE

      set.seed(0)
      if (ppopnw == TRUE) {
        netest_ergm.ego <- netest(nwe, formation = ~edges + degree(1) + nodemix("race"),
                                  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                                  set.control.ergm.ego = control.ergm.ego(ppopsize = nw))
      } else {
        ## assumes node order is preserved in nw -> as.egor(nw) -> pseudopopulation in ergm.ego
        netest_ergm.ego <- netest(nwe, formation = ~edges + degree(1) + nodemix("race"),
                                  coef.diss = dissolution_coefs(~offset(edges), 10, 0))
      }

      set.seed(0)
      netest_ergm <- netest(nw, formation = ~edges + degree(1) + nodemix("race"),
                            target.stats = netest_ergm.ego$target.stats,
                            coef.diss = dissolution_coefs(~offset(edges), 10, 0))

      if (trim == TRUE) {
        netest_ergm <- trim_netest(netest_ergm)
        netest_ergm.ego <- trim_netest(netest_ergm.ego)
      }

      expect_identical(netest_ergm$coef.form.crude, netest_ergm.ego$coef.form.crude)

      set.seed(0)
      dx_ergm <- netdx(netest_ergm, nsteps = 5, dynamic = FALSE, verbose = FALSE)
      set.seed(0)
      dx_ergm.ego <- netdx(netest_ergm.ego, nsteps = 5, dynamic = FALSE, verbose = FALSE)

      expect_identical(dx_ergm$stats, dx_ergm.ego$stats)
      expect_identical(dx_ergm$stats.table.formation, dx_ergm.ego$stats.table.formation)

      set.seed(0)
      dxd_ergm <- netdx(netest_ergm, nsteps = 5, dynamic = TRUE, verbose = FALSE)
      set.seed(0)
      dxd_ergm.ego <- netdx(netest_ergm.ego, nsteps = 5, dynamic = TRUE, verbose = FALSE)

      expect_identical(dxd_ergm$stats, dxd_ergm.ego$stats)
      expect_identical(dxd_ergm$stats.table.formation, dxd_ergm.ego$stats.table.formation)
      expect_identical(dxd_ergm$stats.table.dissolution, dxd_ergm.ego$stats.table.dissolution)
      expect_identical(dxd_ergm$pages, dxd_ergm.ego$pages)
      expect_identical(dxd_ergm$pages_imptd, dxd_ergm.ego$pages_imptd)

      param <- param.net(inf.prob = 0.3, act.rate = 0.5)
      init <- init.net(i.num = 10)
      control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)

      set.seed(0)
      sim_ergm <- netsim(netest_ergm, param, init, control)
      set.seed(0)
      sim_ergm.ego <- netsim(netest_ergm.ego, param, init, control)

      expect_identical(sim_ergm$stats, sim_ergm.ego$stats)
      expect_identical(sim_ergm$epi, sim_ergm.ego$epi)
    }
  }
})
