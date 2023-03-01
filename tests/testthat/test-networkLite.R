context("networkLite (All SOC)")

test_that("network and networkLite work equally in netest, netdx, and netsim", {

  net_size <- 100
  bip_size <- 40

  ffdir <- ~odegree(1) + idegree(1)
  ffundir <- ~degree(1)

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      for (duration in c(1, 10)) {

        nw <- network.initialize(n = 100, directed = directed, bipartite = bipartite)
        nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

        set.seed(0)
        est <- netest(nw, formation = ~edges + nodematch("race"),
                      target.stats = c(50, 20),
                      coef.diss = dissolution_coefs(~offset(edges), c(duration)),
                      verbose = FALSE
        )
        dxs <- netdx(est, nsims = 20, verbose = FALSE,
                       dynamic = FALSE, nwstats.formula = if(directed) ffdir else ffundir)

        if (duration > 1) {
          dxd <- netdx(est, nsims = 2, nsteps = 10, verbose = FALSE,
                         dynamic = TRUE)
        }

        param <- param.net(inf.prob = 0.3, act.rate = 0.5)
        init <- init.net(i.num = 10)
        control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
        sim <- netsim(est, param, init, control)

        nwL <- as.networkLite(nw)
        set.seed(0)
        estL <- netest(nwL, formation = ~edges + nodematch("race"),
                       target.stats = c(50, 20),
                       coef.diss = dissolution_coefs(~offset(edges), c(duration)),
                       verbose = FALSE
        )
        dxsL <- netdx(estL, nsims = 20, verbose = FALSE,
                       dynamic = FALSE, nwstats.formula = if(directed) ffdir else ffundir)

        if (duration > 1) {
          dxdL <- netdx(estL, nsims = 2, nsteps = 10, verbose = FALSE,
                         dynamic = TRUE)
        }

        simL <- netsim(estL, param, init, control)

        # convert networks to networkLites
        dxs$nw <- as.networkLite(dxs$nw)
        if (duration > 1) {
          dxd$nw <- as.networkLite(dxd$nw)
        }

        # the rest should be equal, including coefs, stats, etc.
        expect_equal(trim_netest(est), trim_netest(estL))
        expect_equal(dxs, dxsL)
        if (duration > 1) {
          expect_equal(dxd, dxdL)
        }
        expect_equal(sim, simL)
      }
    }
  }
})
