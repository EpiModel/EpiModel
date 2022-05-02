
## test_that("network and networkLite behave identically in ergm and gof", {
##   skip_on_cran()
##   options(ergm.loglik.warn_dyads=FALSE)
## 
##   net_size <- 100
##   bip_size <- 40
## 
##   for(directed in list(FALSE, TRUE)) {
##     for(bipartite in list(FALSE, bip_size)) {
##       if(directed && bipartite) {
##         next
##       }
## 
##       set.seed(0)
##       nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
##       nw %v% "a" <- rep(letters[1:5], length.out = net_size)
##       nw %v% "b" <- runif(net_size)
##       nw %v% "sex" <- rep(c("M","F"), length.out=net_size)
##       
##       nwL <- as.networkLite(nw)
##       
##       di_constraints <- ~blocks(~sex, levels2=diag(TRUE,2))
##       dd_constraints <- ~bd(maxout=2) + blocks(~sex, levels2=diag(TRUE,2))
##       dm_constraints <- ~bd(maxout=2, minout = 0) + blocks(~sex, levels2=diag(TRUE,2))
##       
##       target_stats <- c(750, 300, 315, 285, 295, 1250)/10
##       
##       set.seed(0)
##       nw_di_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = di_constraints, eval.loglik = FALSE)
##       set.seed(0)
##       nwL_di_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = di_constraints, eval.loglik = FALSE)
##       expect_equal(coef(nw_di_ergm), coef(nwL_di_ergm))
##       
##       set.seed(0)
##       nw_di_gof <- gof(nw_di_ergm)
##       set.seed(0)
##       nwL_di_gof <- gof(nwL_di_ergm)
##       expect_equal(nw_di_gof, nwL_di_gof)
##       
##       set.seed(0)
##       nw_di_predict <- predict(nw_di_ergm)
##       set.seed(0)
##       nwL_di_predict <- predict(nwL_di_ergm)
##       expect_identical(nw_di_predict, nwL_di_predict)
## 
##       set.seed(0)
##       nw_dd_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dd_constraints, control = list(init.method="MPLE"), eval.loglik = FALSE)
##       set.seed(0)
##       nwL_dd_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dd_constraints, control = list(init.method="MPLE"), eval.loglik = FALSE)
##       expect_equal(coef(nw_dd_ergm), coef(nwL_dd_ergm))
## 
##       set.seed(0)
##       nw_dd_gof <- gof(nw_dd_ergm)
##       set.seed(0)
##       nwL_dd_gof <- gof(nwL_dd_ergm)
##       expect_equal(nw_dd_gof, nwL_dd_gof)
##       
##       set.seed(0)
##       nw_dd_predict <- predict(nw_dd_ergm)
##       set.seed(0)
##       nwL_dd_predict <- predict(nwL_dd_ergm)
##       expect_identical(nw_dd_predict, nwL_dd_predict)
## 
##       set.seed(0)
##       nw_dm_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dm_constraints, eval.loglik = FALSE)
##       set.seed(0)
##       nwL_dm_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dm_constraints, eval.loglik = FALSE)
##       expect_equal(coef(nw_dm_ergm), coef(nwL_dm_ergm))
## 
##       set.seed(0)
##       nw_dm_gof <- gof(nw_dm_ergm)
##       set.seed(0)
##       nwL_dm_gof <- gof(nwL_dm_ergm)
##       expect_equal(nw_dm_gof, nwL_dm_gof)
##       
##       set.seed(0)
##       nw_dm_predict <- predict(nw_dm_ergm)
##       set.seed(0)
##       nwL_dm_predict <- predict(nwL_dm_ergm)
##       expect_identical(nw_dm_predict, nwL_dm_predict)
## 
##       ## simpler dyad-independent case where we can hit targets exactly
##       set.seed(0)
##       nw_mple_ergm <- ergm(nw ~ edges + nodefactor("a"), target.stats = as.integer(target_stats[-length(target_stats)]), constraints = di_constraints)
##       set.seed(0)
##       nwL_mple_ergm <- ergm(nwL ~ edges + nodefactor("a"), target.stats = as.integer(target_stats[-length(target_stats)]), constraints = di_constraints)
##       expect_equal(coef(nw_mple_ergm), coef(nwL_mple_ergm))
## 
##       set.seed(0)
##       nw_mple_gof <- gof(nw_mple_ergm)
##       set.seed(0)
##       nwL_mple_gof <- gof(nwL_mple_ergm)
##       expect_equal(nw_mple_gof, nwL_mple_gof)
##       
##       set.seed(0)
##       nw_mple_predict <- predict(nw_mple_ergm)
##       set.seed(0)
##       nwL_mple_predict <- predict(nwL_mple_ergm)
##       expect_identical(nw_mple_predict, nwL_mple_predict)
##       
##     }
##   }
## })

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

test_that("direct conversion between network and networkLite functions as expected", {
  net_size <- 100
  bip_size <- 40

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      for(last.mode in list(FALSE, TRUE)) {
        for(delete in list(FALSE, TRUE)) {
          set.seed(0)
          nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
          nw %v% "b" <- runif(net_size)
          nw <- san(nw ~ edges, target.stats = c(100))
          nw %e% "eattr" <- runif(network.edgecount(nw))
          nw %n% "nattr" <- "attr"
          nainds <- sample(valid.eids(nw), as.integer(network.size(nw)/2), FALSE)
          set.edge.attribute(nw, "na", TRUE, nainds)
          add.vertices(nw, 9, vattr = rep(list(list(na = FALSE, vertex.names = NA_integer_, b = NA_real_)), 9), last.mode = last.mode)
          if(delete) {
            el <- as.edgelist(nw, attrname = "na", na.rm = FALSE)
            w1 <- sample(which(as.logical(el[,3])))[1:5]
            w2 <- sample(which(!as.logical(el[,3])))[1:7]
            delete.edges(nw, unlist(get.dyads.eids(nw, el[w1,1], el[w1,2], na.omit = FALSE)))
            delete.edges(nw, unlist(get.dyads.eids(nw, el[w2,1], el[w2,2], na.omit = FALSE)))
            vd <- sample(seq_len(net_size), 10, FALSE)
            delete.vertices(nw, vd)
          }
          # add.vertices and delete.vertices convert network size to integer....
          nw %n% "n" <- as.numeric(nw %n% "n")
          
          set.seed(0)
          nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
          nwL %v% "b" <- runif(net_size)
          nwL <- san(nwL ~ edges, target.stats = c(100))
          nwL %e% "eattr" <- runif(network.edgecount(nwL))
          nwL %n% "nattr" <- "attr"
          nainds <- sample(valid.eids(nwL), as.integer(network.size(nwL)/2), FALSE)
          set.edge.attribute(nwL, "na", TRUE, nainds)
          add.vertices(nwL, 9, vattr = rep(list(list(na = FALSE, vertex.names = NA_integer_, b = NA_real_)), 9), last.mode = last.mode)
          if(delete) {
            el <- as.edgelist(nwL, attrname = "na", na.rm = FALSE)
            w1 <- sample(which(as.logical(el[,3])))[1:5]
            w2 <- sample(which(!as.logical(el[,3])))[1:7]
            delete.edges(nwL, c(w1,w2))
            vd <- sample(seq_len(net_size), 10, FALSE)
            delete.vertices(nwL, vd)
          }
          
          expect_identical(as.networkLite(nw), nwL)
          expect_identical(as.networkLite(is.na(nw)), is.na(nwL))
          expect_identical(as.networkLite(is.na(is.na(nw))), is.na(is.na(nwL)))
          
          if(delete) {
            expect_identical(as.networkLite(nw), as.networkLite(to_network_networkLite(nwL)))
            expect_identical(as.networkLite(is.na(nw)), as.networkLite(to_network_networkLite(is.na(nwL))))
          } else {
            expect_identical(nw, to_network_networkLite(nwL))
            expect_identical(is.na(nw), to_network_networkLite(is.na(nwL)))
          }
        }
      }
    }
  }
})

## test_that("network and networkLite estimate equally in (EGMME) tergm", {
##   skip_on_cran()
##   net_size <- 50
##   bip_size <- 20
## 
##   for(directed in list(FALSE, TRUE)) {
##     for(bipartite in list(FALSE, bip_size)) {
##       if(directed && bipartite) {
##         next
##       }
## 
##       set.seed(0)
##       nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
##       nw %v% "a" <- rep(letters[1:5], length.out = net_size)
##       nw %v% "b" <- runif(net_size)
##       
##       nwL <- as.networkLite(nw)
##       
##       set.seed(0)
##       tergm_nw <- tergm(nw ~ Form(~edges) + Diss(~edges), targets = ~edges + mean.age, target.stats = c(30, 5), estimate = "EGMME")
##       set.seed(0)
##       tergm_nwL <- tergm(nwL ~ Form(~edges) + Diss(~edges), targets = ~edges + mean.age, target.stats = c(30, 5), estimate = "EGMME")
##       expect_equal(coef(tergm_nw), coef(tergm_nwL))
##     }
##   }
## })

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

      nw <- network_initialize(n = 100, directed = directed, bipartite = bipartite)
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
      
      # drop some things that aren't equal (due to network != networkLite)
      est$fit <- NULL
      estL$fit <- NULL
      
      dxs$nw <- NULL
      dxsL$nw <- NULL
      
      dxd$nw <- NULL
      dxdL$nw <- NULL
      
      # the rest should be equal, including coefs, stats, etc.
      expect_equal(est, estL)
      expect_equal(dxs, dxsL)
      expect_equal(dxd, dxdL)
      expect_equal(sim, simL)
    }
  }
})

test_that("network and networkLite produce identical matrices, edgelists, and tibbles", {
  net_size <- 100
  bip_size <- 40
  edges_target <- net_size
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw <- san(nw ~ edges, target.stats = c(edges_target))
      nw %e% "eattr" <- runif(network.edgecount(nw))
      nainds <- sample(valid.eids(nw), as.integer(length(valid.eids(nw))/2), FALSE)
      set.edge.attribute(nw, "na", TRUE, nainds)

      set.seed(0)
      nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
      nwL <- san(nwL ~ edges, target.stats = c(edges_target))
      nwL %e% "eattr" <- runif(network.edgecount(nwL))
      nainds <- sample(valid.eids(nwL), as.integer(length(valid.eids(nwL))/2), FALSE)
      set.edge.attribute(nwL, "na", TRUE, nainds)
      
      for(attrname in list(NULL, "eattr", "na")) {
        for(na.rm in list(FALSE, TRUE)) {
          for(matrix.type in c("adjacency", "incidence", "edgelist")) {
            expect_identical(as.matrix(nw, matrix.type = matrix.type, attrname = attrname, na.rm = na.rm),
                             as.matrix(nwL, matrix.type = matrix.type, attrname = attrname, na.rm = na.rm))
          }
          expect_identical(as.edgelist(nw, attrname = attrname, na.rm = na.rm),
                           as.edgelist(nwL, attrname = attrname, na.rm = na.rm))

          expect_identical(as.edgelist(nw, attrname = attrname, na.rm = na.rm, output = "tibble"),
                           as.edgelist(nwL, attrname = attrname, na.rm = na.rm, output = "tibble"))

          expect_identical(tibble::as_tibble(nw, attrname = attrname, na.rm = na.rm),
                           tibble::as_tibble(nwL, attrname = attrname, na.rm = na.rm))
        }
      }
    }
  }
})

## test_that("network and networkLite fit and simulate equal missing-data ergms", {
##   net_size <- 50
##   bip_size <- 20
##   
##   for(directed in list(FALSE, TRUE)) {
##     for(bipartite in list(FALSE, bip_size)) {
##       if(directed && bipartite) {
##         next
##       }
##       
##       set.seed(0)
##       nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
##       nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL)/10)
##       nwL %v% "age" <- runif(net_size)
##       na <- sample(c(FALSE,TRUE),network.edgecount(nwL),TRUE)
##       
##       set.seed(0)
##       eL <- ergm(nwL ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
##       nwL %e% "na" <- na
##       set.seed(0)
##       eLna <- ergm(nwL ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
##       eL2 <- simulate(eLna)
##       
##       set.seed(0)
##       nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
##       nw <- san(nw ~ edges, target.stats = network.dyadcount(nw)/10)
##       nw %v% "age" <- runif(net_size)
##       na <- sample(c(FALSE,TRUE),network.edgecount(nw),TRUE)
## 
##       set.seed(0)
##       e <- ergm(nw ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
##       nw %e% "na" <- na
##       set.seed(0)
##       ena <- ergm(nw ~ absdiff("age"), control = list(MCMLE.effectiveSize = NULL))
##       e2 <- simulate(ena)
##       
##       expect_equal(coef(e), coef(eL))
##       expect_equal(coef(ena), coef(eLna))
##       expect_equal(as.edgelist(e2), as.edgelist(eL2))
##       expect_equal(as.edgelist(e2, attrname = "na"), as.edgelist(eL2, attrname = "na"))
##     }
##   }
## })

## test_that("network and networkLite fit and simulate equal valued ergms", {
##   net_size <- 50
##   bip_size <- 20
##   
##   for(directed in list(FALSE, TRUE)) {
##     for(bipartite in list(FALSE, bip_size)) {
##       if(directed && bipartite) {
##         next
##       }
##       
##       set.seed(0)
##       nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
##       nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL))
##       nwL %v% "age" <- runif(net_size)
##       nwL %e% "w" <- runif(network.edgecount(nwL))
##       eL <- ergm(nwL ~ absdiff("age"), response = "w", reference = ~Unif(0,1), control = list(MCMLE.effectiveSize = NULL))
##       eL2 <- simulate(eL)
##       
##       set.seed(0)
##       nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
##       nw <- san(nw ~ edges, target.stats = network.dyadcount(nw))
##       nw %v% "age" <- runif(net_size)
##       nw %e% "w" <- runif(network.edgecount(nw))
##       e <- ergm(nw ~ absdiff("age"), response = "w", reference = ~Unif(0,1), control = list(MCMLE.effectiveSize = NULL))
##       e2 <- simulate(e)
##       
##       expect_equal(coef(e), coef(eL))
##       expect_equal(as.edgelist(e2, attrname = "w"), as.edgelist(eL2, attrname = "w"))
##     }
##   }
## })

test_that("network and networkLite `[<-` and add.edges produce consistent edgelists", {
  net_size <- 100
  bip_size <- 40
  edges_target <- net_size
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      b1 <- if(bipartite) bip_size else net_size
      b2 <- if(bipartite) net_size - bip_size else net_size
      
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)

      nwa <- nw
      nwLa <- nwL

      nw0 <- nw

      rv <- san(nw0 ~ edges, target.stats = c(edges_target))
      m <- as.matrix(rv, matrix.type = "adjacency")
      el <- as.matrix(rv, matrix.type = "edgelist")
      
      nw[,] <- m
      nwL[,] <- m
      nwa <- add.edges(nwa, el[,1], el[,2])
      nwLa <- add.edges(nwLa, el[,1], el[,2])
      
      expect_equal(as.edgelist(nw), as.edgelist(nwL))
      expect_equal(as.edgelist(nwa), as.edgelist(nwLa))

      rv2 <- san(nw0 ~ edges, target.stats = c(edges_target))
      m2 <- as.matrix(rv2, matrix.type = "adjacency")
      el2 <- as.matrix(rv2 - rv, matrix.type = "edgelist")

      nw[,] <- m2
      nwL[,] <- m2
      nwa <- add.edges(nwa, el2[,1], el2[,2])
      nwLa <- add.edges(nwLa, el2[,1], el2[,2])
      
      expect_equal(as.edgelist(nw), as.edgelist(nwL))
      expect_equal(as.edgelist(nwa), as.edgelist(nwLa))

      m <- matrix(runif(b1*b2), b1, b2)
      if(!directed && !bipartite) {
        m <- m + t(m)
      }
      
      m[m < 0.5] <- 0
      nwm <- network(m > 0, matrix.type = if(bipartite) "bipartite" else "adjacency", directed = directed, bipartite = bipartite)
      nwm[,,names.eval="w",add.edges=FALSE] <- m
      
      nwmd <- nwm - (rv + rv2)
      nwmd[,,names.eval="w",add.edges=FALSE] <- m
      elm <- as.edgelist(nwmd, attrname="w")
      
      nw[,,names.eval="w",add.edges=FALSE] <- m
      nwL[,,names.eval="w",add.edges=FALSE] <- m      
      
      expect_equal(as.edgelist(nw, attrname = "w"), as.edgelist(nwL, attrname = "w"))
      
      nw[,,names.eval="w",add.edges=TRUE] <- m
      nwL[,,names.eval="w",add.edges=TRUE] <- m
      nwa <- add.edges(nwa, elm[,1], elm[,2], names.eval="w", vals.eval=elm[,3])
      nwLa <- add.edges(nwLa, elm[,1], elm[,2], names.eval="w", vals.eval=elm[,3])

      expect_equal(as.edgelist(nw, attrname = "w"), as.edgelist(nwL, attrname = "w"))
      expect_equal(as.edgelist(nwa, attrname = "w"), as.edgelist(nwLa, attrname = "w"))      
      
      nw[,] <- FALSE
      nwL[,] <- FALSE
      nwa[,] <- FALSE
      nwLa[,] <- FALSE

      elm2 <- as.edgelist(nwm, attrname="w")
      
      nw[,,names.eval="w",add.edges=TRUE] <- m
      nwL[,,names.eval="w",add.edges=TRUE] <- m

      nwa <- add.edges(nwa, elm2[,1], elm2[,2], names.eval="w", vals.eval=elm2[,3])
      nwLa <- add.edges(nwLa, elm2[,1], elm2[,2], names.eval="w", vals.eval=elm2[,3])
      
      expect_equal(as.edgelist(nw, attrname = "w"), as.edgelist(nwL, attrname = "w"))
      expect_equal(as.edgelist(nwa, attrname = "w"), as.edgelist(nwLa, attrname = "w"))      
    }
  }  
})


test_that("network and networkLite `+` and `-` produce consistent results", {
  net_size <- 100
  bip_size <- 40
  edges_target <- 10*net_size
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      
      nw1 <- san(nw ~ edges, target.stats = c(edges_target))
      nw2 <- san(nw ~ edges, target.stats = c(edges_target))

      nwL1 <- as.networkLite(nw1)
      nwL2 <- as.networkLite(nw2)
      
      expect_identical(as.edgelist(nw1 + nw2), as.edgelist(nwL1 + nwL2))
      expect_identical(as.edgelist(nw1 - nw2), as.edgelist(nwL1 - nwL2))
    }
  }
})

test_that("network to networkLite conversion handles deleted edges with attributes appropriately", {
  net_size <- 10
  bip_size <- 4
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw[1,5] <- 1
      nw[4,7] <- 1
      nw[3,9] <- 1
      nw[2,6] <- 1
      nw[1,10] <- 1
      eattr <- runif(5)
      nw %e% "eattr" <- eattr
      delete.edges(nw, c(2,4))
      
      nwL <- as.networkLite(nw)
      expect_identical(nw %e% "eattr", eattr[c(1,3,5)])
      expect_identical(nwL %e% "eattr", eattr[c(1,5,3)])
      expect_identical(as.edgelist(nw), as.edgelist(nwL))
    }
  }
})

test_that("network and networkLite behave equivalently for basic access and mutation", {
  net_size <- 10
  bip_size <- 4
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw[1,5] <- 1
      nw[4,7] <- 1
      nw[3,9] <- 1
      nw[2,6] <- 1
      nw[1,10] <- 1
      eattr1 <- runif(3)
      e1 <- c(2,1,4)
      eattr2 <- "a"
      e2 <- c(2,4)
      eattr3 <- runif(5)
      
      set.edge.attribute(nw, "eattr1", eattr1, e1)
      set.edge.attribute(nw, "eattr2", eattr2, e2)
      nw %e% "eattr3" <- eattr3

      vattr1 <- sample(c("a","b"), 7, TRUE)
      v1 <- c(5,9,4,6,8,1,2)
      vattr2 <- FALSE
      v2 <- c(8)
      vattr3 <- runif(10)
            
      set.vertex.attribute(nw, "vattr1", vattr1, v1)
      set.vertex.attribute(nw, "vattr2", vattr2, v2)
      nw %v% "vattr3" <- vattr3
      
      nwL <- networkLite(as.edgelist(nw))
      
      eo <- c(1,5,4,3,2)
      
      set.edge.attribute(nwL, "eattr1", eattr1, eo[e1])
      set.edge.attribute(nwL, "eattr2", eattr2, eo[e2])
      nwL %e% "eattr3" <- eattr3[eo]

      set.vertex.attribute(nwL, "vattr1", vattr1, v1)
      set.vertex.attribute(nwL, "vattr2", vattr2, v2)
      nwL %v% "vattr3" <- vattr3
      
      expect_identical(as.edgelist(nw, attrname = "eattr1"), as.edgelist(nwL, attrname = "eattr1"))
      expect_identical(as.edgelist(nw, attrname = "eattr2"), as.edgelist(nwL, attrname = "eattr2"))
      expect_identical(as.edgelist(nw, attrname = "eattr3"), as.edgelist(nwL, attrname = "eattr3"))
      
      expect_identical(nw %v% "vattr1", nwL %v% "vattr1")
      expect_identical(nw %v% "vattr2", nwL %v% "vattr2")
      expect_identical(nw %v% "vattr3", nwL %v% "vattr3")
    }
  }
})

test_that("add.vertices and add.edges with irregular attribute arguments behave equivalently for network and networkLite", {
  net_size <- 100
  bip_size <- 40

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      for(last.mode in list(FALSE, TRUE)) {
        
        vnames <- paste0("v", 1:4)
        enames <- paste0("e", 1:4)
        
        set.seed(0)
        nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
        nwe <- san(nw ~ edges, target.stats = c(net_size))
        nw <- san(nw ~ edges, target.stats = c(net_size))
        nw %v% "v1" <- runif(net_size)
        nw %v% "v2" <- runif(net_size)
        nw %e% "e1" <- runif(network.edgecount(nw))
        nw %e% "e2" <- runif(network.edgecount(nw))
        el <- as.edgelist(nwe - nw)
        names.eval <- list()
        vals.eval <- list()
        for(i in seq_len(NROW(el))) {
          en <- which(as.logical(round(runif(length(enames)))))
          if(length(en) > 0) {
            en <- sample(en)
            names.eval[[i]] <- as.list(enames[en])
            vals.eval[[i]] <- as.list(runif(length(en)))
          } else {
            names.eval[[i]] <- list()
            vals.eval[[i]] <- list()
          }
        }
        add.edges(nw, el[,1], el[,2], names.eval = names.eval, vals.eval = vals.eval)
        
        vta <- 50
        vattr <- list()
        for(i in seq_len(vta)) {
          vn <- which(as.logical(round(runif(length(vnames)))))
          if(length(vn) > 0) {
            vn <- sample(vn)
            vattr[[i]] <- as.list(runif(length(vn)))
            names(vattr[[i]]) <- vnames[vn]
          } else {
            vattr[[i]] <- list()
          }
        }
        add.vertices(nw, vta, vattr = vattr, last.mode = last.mode)
        
        set.seed(0)
        nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
        nwLe <- san(nwL ~ edges, target.stats = c(net_size))
        nwL <- san(nwL ~ edges, target.stats = c(net_size))
        nwL %v% "v1" <- runif(net_size)
        nwL %v% "v2" <- runif(net_size)
        nwL %e% "e1" <- runif(network.edgecount(nwL))
        nwL %e% "e2" <- runif(network.edgecount(nwL))
        el <- as.edgelist(nwLe - nwL)

        names.eval <- list()
        vals.eval <- list()
        for(i in seq_len(NROW(el))) {
          en <- which(as.logical(round(runif(length(enames)))))
          if(length(en) > 0) {
            en <- sample(en)
            names.eval[[i]] <- as.list(enames[en])
            vals.eval[[i]] <- as.list(runif(length(en)))
          } else {
            names.eval[[i]] <- list()
            vals.eval[[i]] <- list()
          }
        }
        add.edges(nwL, el[,1], el[,2], names.eval = names.eval, vals.eval = vals.eval)
        
        vta <- 50
        vattr <- list()
        for(i in seq_len(vta)) {
          vn <- which(as.logical(round(runif(length(vnames)))))
          if(length(vn) > 0) {
            vn <- sample(vn)
            vattr[[i]] <- as.list(runif(length(vn)))
            names(vattr[[i]]) <- vnames[vn]
          } else {
            vattr[[i]] <- list()
          }
        }
        add.vertices(nwL, vta, vattr = vattr, last.mode = last.mode)

        for(en in setdiff(list.edge.attributes(nw), "na")) {
          ev <- get.edge.attribute(nwL, en)
          delete.edge.attribute(nwL, en)
          set.edge.attribute(nwL, en, ev)
        }
        
        for(vn in setdiff(list.vertex.attributes(nw), c("na", "vertex.names"))) {
          vv <- get.vertex.attribute(nwL, vn)
          delete.vertex.attribute(nwL, vn)
          set.vertex.attribute(nwL, vn, vv)
        }

        expect_equal(as.networkLite(nw), nwL)
        expect_equal(as.networkLite(nw), as.networkLite(to_network_networkLite(nwL)))
      }
    }
  }
})

test_that("attribute setting and deleting behave equivalently for network and networkLite", {
  net_size <- 10
  bip_size <- 4

  enames <- paste0("e", 1:10)
  vnames <- paste0("v", 1:10)
  nnames <- paste0("n", 1:10)
  niter <- 100

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }
      
      set.seed(0)
      nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw <- san(nw ~ edges, target.stats = c(net_size))
      for(i in seq_len(niter)) {
        en <- sample(enames, 1)
        vn <- sample(vnames, 1)
        nn <- sample(nnames, 1)
        
        if(en %in% list.edge.attributes(nw)) {
          delete.edge.attribute(nw, en)
        } else {
          set.edge.attribute(nw, en, runif(network.edgecount(nw)))
        }
        if(vn %in% list.vertex.attributes(nw)) {
          delete.vertex.attribute(nw, vn)
        } else {
          set.vertex.attribute(nw, vn, runif(network.size(nw)))
        }
        if(nn %in% list.network.attributes(nw)) {
          delete.network.attribute(nw, nn)
        } else {
          set.network.attribute(nw, nn, runif(1))
        }
      }
      
      set.seed(0)
      nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
      nwL <- san(nwL ~ edges, target.stats = c(net_size))
      for(i in seq_len(niter)) {
        en <- sample(enames, 1)
        vn <- sample(vnames, 1)
        nn <- sample(nnames, 1)
        
        if(en %in% list.edge.attributes(nwL)) {
          delete.edge.attribute(nwL, en)
        } else {
          set.edge.attribute(nwL, en, runif(network.edgecount(nwL)))
        }
        if(vn %in% list.vertex.attributes(nwL)) {
          delete.vertex.attribute(nwL, vn)
        } else {
          set.vertex.attribute(nwL, vn, runif(network.size(nwL)))
        }
        if(nn %in% list.network.attributes(nwL)) {
          delete.network.attribute(nwL, nn)
        } else {
          set.network.attribute(nwL, nn, runif(1))
        }
      }

      ## re-order everything for these comparisons...
      for(en in setdiff(list.edge.attributes(nw), "na")) {
        ev <- get.edge.attribute(nwL, en)
        delete.edge.attribute(nwL, en)
        set.edge.attribute(nwL, en, ev)
      }
      
      for(vn in setdiff(list.vertex.attributes(nw), c("na", "vertex.names"))) {
        vv <- get.vertex.attribute(nwL, vn)
        delete.vertex.attribute(nwL, vn)
        set.vertex.attribute(nwL, vn, vv)
      }

      for(nn in setdiff(list.network.attributes(nw), c("n", "directed", "bipartite", "loops", "hyper", "multiple", "mnext"))) {
        nv <- get.network.attribute(nwL, nn)
        delete.network.attribute(nwL, nn)
        set.network.attribute(nwL, nn, nv)      
      }

      for(en in setdiff(list.edge.attributes(nwL), "na")) {
        ev <- get.edge.attribute(nw, en)
        delete.edge.attribute(nw, en)
        set.edge.attribute(nw, en, ev)
      }
      
      for(vn in setdiff(list.vertex.attributes(nwL), c("na", "vertex.names"))) {
        vv <- get.vertex.attribute(nw, vn)
        delete.vertex.attribute(nw, vn)
        set.vertex.attribute(nw, vn, vv)
      }

      for(nn in setdiff(list.network.attributes(nwL), c("n", "directed", "bipartite", "loops", "hyper", "multiple", "mnext"))) {
        nv <- get.network.attribute(nw, nn)
        delete.network.attribute(nw, nn)
        set.network.attribute(nw, nn, nv)      
      }
      
      expect_identical(nw, to_network_networkLite(nwL))
      expect_identical(as.networkLite(nw), nwL)
    }
  }
})
