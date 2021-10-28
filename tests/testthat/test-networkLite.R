test_that("network and networkLite estimate identically in ergm", {
  require(ergm)

  options(ergm.loglik.warn_dyads=FALSE)

  set.seed(0)
  nw <- network.initialize(1000, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 200)
  nw %v% "b" <- runif(1000)
  nw %v% "sex" <- rep(c("M","F"), length.out=1000)

  nwL <- as.networkLite(nw)
  
  di_constraints <- ~blocks(~sex, levels2=diag(TRUE,2))
  dd_constraints <- ~bd(maxout=2) + blocks(~sex, levels2=diag(TRUE,2))
  dm_constraints <- ~bd(maxout=2, minout = 0) + blocks(~sex, levels2=diag(TRUE,2))
  
  target_stats <- c(851.0370, 375.2088, 384.6334, 357.3602, 250.4054, 1468.3650)

  set.seed(0)
  nw_di_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = di_constraints, eval.loglik = FALSE)
  set.seed(0)
  nwL_di_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = di_constraints, eval.loglik = FALSE)
  expect_equal(coef(nw_di_ergm), coef(nwL_di_ergm))

  set.seed(0)
  nw_dd_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dd_constraints, control = list(init.method="MPLE"), eval.loglik = FALSE)
  set.seed(0)
  nwL_dd_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dd_constraints, control = list(init.method="MPLE"), eval.loglik = FALSE)
  expect_equal(coef(nw_dd_ergm), coef(nwL_dd_ergm))

  set.seed(0)
  nw_dm_ergm <- ergm(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dm_constraints, eval.loglik = FALSE)
  set.seed(0)
  nwL_dm_ergm <- ergm(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = target_stats, constraints = dm_constraints, eval.loglik = FALSE)
  expect_equal(coef(nw_dm_ergm), coef(nwL_dm_ergm))

  ## simpler dyad-independent case where we can hit targets exactly  
  set.seed(0)
  nw_mple_ergm <- ergm(nw ~ edges + nodefactor("a"), target.stats = as.integer(target_stats[-length(target_stats)]), constraints = di_constraints)
  set.seed(0)
  nwL_mple_ergm <- ergm(nwL ~ edges + nodefactor("a"), target.stats = as.integer(target_stats[-length(target_stats)]), constraints = di_constraints)
  expect_equal(coef(nw_mple_ergm), coef(nwL_mple_ergm))
})


test_that("network and networkLite simulate identically in ergm", {
  require(ergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
  
  coef <- c(-4, 1, 1.5, 0.5, -1, 0.5)
  
  set.seed(0)
  nw_1 <- simulate(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  set.seed(0)
  nwL_1 <- simulate(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)  
  expect_identical(summary(nw_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)), 
                   summary(nwL_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)))
  
  set.seed(0)
  nw_2 <- simulate(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  set.seed(0)
  nwL_2 <- simulate(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)  
  expect_identical(summary(nw_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)), 
                   summary(nwL_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)))
})

test_that("network and networkLite simulate identically in san", {
  require(ergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
    
  set.seed(0)
  nw_1 <- san(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
  set.seed(0)
  nwL_1 <- san(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)
  expect_identical(summary(nw_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)), 
                   summary(nwL_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)))
  
  set.seed(0)
  nw_2 <- san(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
  set.seed(0)
  nwL_2 <- san(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)
  expect_identical(summary(nw_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)), 
                   summary(nwL_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + gwnsp(0.3, fixed=TRUE)))
})

test_that("network and networkLite simulate identically in tergm", {
  require(tergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
  
  coef <- c(-4, 1, 1.5, 0.5, -1, 0.5, 3)
  
  set.seed(0)
  nw_1 <- simulate(nw ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_1 <- simulate(nwL ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)
  expect_identical(nw_1 %n% "lasttoggle", nwL_1 %n% "lasttoggle")
  expect_identical(nw_1 %n% "time", nwL_1 %n% "time")
  expect_identical(summary(nw_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)),
                   summary(nwL_1 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)))
  
  set.seed(0)
  nw_2 <- simulate(nw_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_2 <- simulate(nwL_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)  
  expect_identical(nw_2 %n% "lasttoggle", nwL_2 %n% "lasttoggle")
  expect_identical(nw_2 %n% "time", nwL_2 %n% "time")
  expect_identical(summary(nw_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)),
                   summary(nwL_2 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)))

  set.seed(0)
  nw_3 <- simulate(nw_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_3 <- simulate(nwL_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Persist(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_3)), unclass(as.edgelist(nwL_3)), check.attributes = FALSE)  
  expect_identical(nw_3 %n% "lasttoggle", nwL_3 %n% "lasttoggle")
  expect_identical(nw_3 %n% "time", nwL_3 %n% "time")
  expect_identical(summary(nw_3 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)),
                   summary(nwL_3 ~ nodemix(~a) + absdiff(~b) + concurrent + gwesp + mean.age + edge.ages + nodemix.mean.age(~a) + gwnsp(0.3, fixed=TRUE)))
})
