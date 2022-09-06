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

      nw <- network.initialize(n = 100, directed = directed, bipartite = bipartite)
      nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

      set.seed(0)
      est <- netest(nw, formation = ~edges + nodematch("race"),
                    target.stats = c(50, 20),
                    coef.diss = dissolution_coefs(~offset(edges), c(10)),
                    verbose = FALSE
      )
      dxs <- netdx(est, nsims = 20, verbose = FALSE,
                     dynamic = FALSE, nwstats.formula = if(directed) ffdir else ffundir)

      dxd <- netdx(est, nsims = 2, nsteps = 10, verbose = FALSE,
                     dynamic = TRUE)

      param <- param.net(inf.prob = 0.3, act.rate = 0.5)
      init <- init.net(i.num = 10)
      control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
      sim <- netsim(est, param, init, control)

      nwL <- as.networkLite(nw)
      set.seed(0)
      estL <- netest(nwL, formation = ~edges + nodematch("race"),
                     target.stats = c(50, 20),
                     coef.diss = dissolution_coefs(~offset(edges), c(10)),
                     verbose = FALSE
      )
      dxsL <- netdx(estL, nsims = 20, verbose = FALSE,
                     dynamic = FALSE, nwstats.formula = if(directed) ffdir else ffundir)

      dxdL <- netdx(estL, nsims = 2, nsteps = 10, verbose = FALSE,
                     dynamic = TRUE)

      simL <- netsim(estL, param, init, control)

      # convert networks to networkLites
      dxs$nw <- as.networkLite(dxs$nw)
      dxd$nw <- as.networkLite(dxd$nw)

      # the rest should be equal, including coefs, stats, etc.
      expect_equal(trim_netest(est), trim_netest(estL))
      expect_equal(dxs, dxsL)
      expect_equal(dxd, dxdL)
      expect_equal(sim, simL)
    }
  }
})

test_that("network and networkLite simulate equally in tergm", {
  net_size <- 100
  bip_size <- 40

  ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)
  ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw %v% "a" <- rep(letters[1:5], length.out = net_size)
      nw %v% "b" <- runif(net_size)

      nwL <- as.networkLite(nw)

      coef <- c(-4, 1, 1.5, 0.5, -1, 0.5, 3)

      set.seed(0)
      nw_1 <- simulate(nw ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      set.seed(0)
      nwL_1 <- simulate(nwL ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      expect_is(nwL_1, "networkLite")

      expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
      expect_identical(nw_1 %n% "lasttoggle", nwL_1 %n% "lasttoggle")
      expect_identical(nw_1 %n% "time", nwL_1 %n% "time")
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_1),
                         summary(ffdir, basis = nwL_1))
      } else {
        expect_identical(summary(ffundir, basis = nw_1),
                         summary(ffundir, basis = nwL_1))
      }

      set.seed(0)
      nw_2 <- simulate(nw_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      set.seed(0)
      nwL_2 <- simulate(nwL_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      expect_is(nwL_2, "networkLite")

      expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
      expect_identical(nw_2 %n% "lasttoggle", nwL_2 %n% "lasttoggle")
      expect_identical(nw_2 %n% "time", nwL_2 %n% "time")
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_2),
                         summary(ffdir, basis = nwL_2))
      } else {
        expect_identical(summary(ffundir, basis = nw_2),
                         summary(ffundir, basis = nwL_2))
      }

      set.seed(0)
      nw_3 <- simulate(nw_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      set.seed(0)
      nwL_3 <- simulate(nwL_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
      expect_is(nwL_3, "networkLite")

      expect_equal(as.edgelist(nw_3), as.edgelist(nwL_3))
      expect_identical(nw_3 %n% "lasttoggle", nwL_3 %n% "lasttoggle")
      expect_identical(nw_3 %n% "time", nwL_3 %n% "time")
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_3),
                         summary(ffdir, basis = nwL_3))
      } else {
        expect_identical(summary(ffundir, basis = nw_3),
                         summary(ffundir, basis = nwL_3))
      }

      set.seed(0)
      nw_4 <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE)
      set.seed(0)
      nwL_4 <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE)

      # comparison of networkDynamics
      expect_equal(nw_4, nwL_4)


      ## for completeness, also get stats and changes as output
      set.seed(0)
      s <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE, output = "stats", stats = TRUE, monitor = if(directed) ~edges + idegree(0:10) + odegree(0:10) + mean.age + Form(~odegree(0:2)) else ~edges + degree(0:10) + mean.age + Form(~degree(0:2)))
      set.seed(0)
      sL <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE, output = "stats", stats = TRUE, monitor = if(directed) ~edges + idegree(0:10) + odegree(0:10) + mean.age + Form(~odegree(0:2)) else ~edges + degree(0:10) + mean.age + Form(~degree(0:2)))

      # comparison of stats
      expect_equal(s, sL)

      set.seed(0)
      c <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE, output = "changes")
      set.seed(0)
      cL <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE, output = "changes")

      # comparison of changes
      expect_equal(c, cL)

      # again, without lasttoggle
      nw_3 %n% "lasttoggle" <- NULL
      nwL_3 %n% "lasttoggle" <- NULL

      set.seed(0)
      nw_4 <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE)
      set.seed(0)
      nwL_4 <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, dynamic = TRUE)

      # comparison of networkDynamics
      expect_equal(nw_4, nwL_4)
    }
  }
})

test_that("network and networkLite simulate equally in ergm", {
  net_size <- 100
  bip_size <- 40

  ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) + gwesp + gwnsp(0.3, fixed=TRUE)
  ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw %v% "a" <- rep(letters[1:5], length.out = net_size)
      nw %v% "b" <- runif(net_size)

      nwL <- as.networkLite(nw)

      coef <- c(-4, 1, 1.5, 0.5, -1, 0.5)

      set.seed(0)
      nw_1 <- simulate(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
      set.seed(0)
      nwL_1 <- simulate(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
      expect_is(nwL_1, "networkLite")

      expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_1),
                         summary(ffdir, basis = nwL_1))
      } else {
        expect_identical(summary(ffundir, basis = nw_1),
                         summary(ffundir, basis = nwL_1))
      }

      set.seed(0)
      nw_2 <- simulate(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
      set.seed(0)
      nwL_2 <- simulate(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
      expect_is(nwL_2, "networkLite")

      expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_2),
                         summary(ffdir, basis = nwL_2))
      } else {
        expect_identical(summary(ffundir, basis = nw_2),
                         summary(ffundir, basis = nwL_2))
      }
    }
  }
})

test_that("network and networkLite simulate equally in san", {
  net_size <- 100
  bip_size <- 40

  ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) + gwesp + gwnsp(0.3, fixed=TRUE)
  ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw %v% "a" <- rep(letters[1:5], length.out = net_size)
      nw %v% "b" <- runif(net_size)

      nwL <- as.networkLite(nw)

      set.seed(0)
      nw_1 <- san(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
      set.seed(0)
      nwL_1 <- san(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
      expect_is(nwL_1, "networkLite")

      expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_1),
                         summary(ffdir, basis = nwL_1))
      } else {
        expect_identical(summary(ffundir, basis = nw_1),
                         summary(ffundir, basis = nwL_1))
      }

      set.seed(0)
      nw_2 <- san(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
      set.seed(0)
      nwL_2 <- san(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
      expect_is(nwL_2, "networkLite")

      expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
      if(directed) {
        expect_identical(summary(ffdir, basis = nw_2),
                         summary(ffdir, basis = nwL_2))
      } else {
        expect_identical(summary(ffundir, basis = nw_2),
                         summary(ffundir, basis = nwL_2))
      }
    }
  }
})

test_that("network and networkLite fit and simulate equal missing-data ergms", {
  net_size <- 50
  bip_size <- 20
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      set.seed(0)
      nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
      nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL)/10)
      nwL %v% "age" <- runif(net_size)
      na <- sample(c(FALSE,TRUE),network.edgecount(nwL),TRUE)
      
      set.seed(0)
      eL <- ergm(nwL ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
      expect_is(eL$newnetwork, "networkLite")
      set.edge.attribute(nwL, "na", na)
      set.seed(0)
      eLna <- ergm(nwL ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
      expect_is(eLna$newnetwork, "networkLite")
      eL2 <- simulate(eLna)
      expect_is(eL2, "networkLite")

      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw <- san(nw ~ edges, target.stats = network.dyadcount(nw)/10)
      nw %v% "age" <- runif(net_size)
      na <- sample(c(FALSE,TRUE),network.edgecount(nw),TRUE)

      set.seed(0)
      e <- ergm(nw ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
      set.edge.attribute(nw, "na", na)
      set.seed(0)
      ena <- ergm(nw ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
      e2 <- simulate(ena)
      
      expect_equal(coef(e), coef(eL))
      expect_equal(coef(ena), coef(eLna))
      expect_equal(as.edgelist(e2), as.edgelist(eL2))
      expect_equal(as.edgelist(e2, attrname = "na"), as.edgelist(eL2, attrname = "na"))
    }
  }
})

test_that("network and networkLite fit and simulate equal valued ergms", {
  net_size <- 50
  bip_size <- 20
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      set.seed(0)
      nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
      nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL))
      nwL %v% "age" <- runif(net_size)
      set.edge.attribute(nwL, "w", runif(network.edgecount(nwL)))
      eL <- ergm(nwL ~ absdiff("age"), response = "w", reference = ~Unif(0,1), control = list(MCMLE.effectiveSize = NULL))
      expect_is(eL$newnetwork, "networkLite")
      eL2 <- simulate(eL)
      expect_is(eL2, "networkLite")

      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw <- san(nw ~ edges, target.stats = network.dyadcount(nw))
      nw %v% "age" <- runif(net_size)
      set.edge.attribute(nw, "w", runif(network.edgecount(nw)))
      e <- ergm(nw ~ absdiff("age"), response = "w", reference = ~Unif(0,1), control = list(MCMLE.effectiveSize = NULL))
      e2 <- simulate(e)
      
      expect_equal(coef(e), coef(eL))
      expect_equal(as.edgelist(e2, attrname = "w"), as.edgelist(eL2, attrname = "w"))
    }
  }
})
