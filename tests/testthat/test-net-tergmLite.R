context("Network models with tergmLite")

# tergmLite = TRUE -----

test_that("tergmLite: 1G, Closed", {

  # SI
  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, act.rate = 5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = 1)
  summary(sim, at = 20)


  # SIS
  param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = 1)
  summary(sim, at = 10)


  # SIR
  param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = 1)
  summary(sim, at = 10)

})

test_that("tergmLite: 2G, Closed", {

  # SI
  num1 <- num2 <- 500
  nw <- network_initialize(n = num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2)
  init <- init.net(i.num = 20, i.num.g2 = 20)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim)
  summary(sim, at = 20)

  # SIS
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 20, i.num.g2 = 20)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim)
  summary(sim, at = 10)

  # SIR
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, mean.col = c(1, 2, 3, 1, 2, 3))
  summary(sim, at = 10)

})


test_that("tergmLite: 1G, Open", {

  # SI
  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10,
                                 d.rate = 0.005)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.4, act.rate = 1,
                     a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE,
                         verbose = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, y = "si.flow", ylim = c(0, 20))
  summary(sim, at = 10)

  # SIS
  param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                     a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, y = c("si.flow", "is.flow"), legend = TRUE)
  summary(sim, at = 10)

  # SIR
  param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                     a.rate = 0.005, di.rate = 0.005, ds.rate = 0.005,
                     dr.rate = 0.005)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1)
  plot(sim, y = "num", ylim = c(800, 1200))

})


test_that("tergmLite: 2G, Open", {

  # SI
  num1 <- num2 <- 500
  nw <- network_initialize(n = num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                                 d.rate = 0.005)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                     act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))

  # SIS
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                     act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))

  # SIR
  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                     act.rate = 5, a.rate = 0.005, a.rate.g2 = 0.005,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                     dr.rate = 0.005, dr.rate.g2 = 0.005,
                     rec.rate = 0.005, rec.rate.g2 = 0.005)
  init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE)

  sim <- netsim(est, param, init, control)
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))

})


test_that("Models with duration = 1", {

  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, act.rate = 5)
  init <- init.net(i.num = 10)
  
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = FALSE, resimulate.network = TRUE)
  sim <- netsim(est, param, init, control)

  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = TRUE, resimulate.network = TRUE)
  sim <- netsim(est, param, init, control)

})

test_that("edgecov works with tergmLite", {
  nw <- network_initialize(n = 20)
  nw %v% "race" <- rep(c("a", "b"), length.out = 20)

  # create matrix for edgecov term
  ec <- matrix(runif(400), 20, 20)
  ec <- (ec + t(ec))/3
  nw %n% "ec" <- ec

  fit <- netest(nw, 
                formation = ~edges + nodefactor("race") + edgecov("ec") , 
                target.stats = c(30, 15, 10), 
                coef.diss = dissolution_coefs(dissolution = ~offset(edges), duration = 1), edapprox = TRUE)

  coefsign <- sign(fit$coef.form.crude[3])

  resim_fun <- function(dat, at) {
    n <- attr(dat$el[[1]], "n")
    m <- matrix(if(at == 3) +Inf else -Inf, n, n)*coefsign
    attr(dat$el[[1]], "ec") <- m
    resim_nets(dat, at)
  }

  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 5, r.num = 0)
  control <- control.net(type = NULL, nsteps = 4, nsims = 1, ncores = 1,
                         resim_nets.FUN = resim_fun, tergmLite = TRUE,
                         resimulate.network = TRUE, nwstats.formula = ~edges,
                         set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5))
  sim <- netsim(fit, param, init, control)

  expect_equal(unname(sim$stats$nwstats$sim1[[1]][2,1,drop=TRUE]), 0)
  expect_equal(unname(sim$stats$nwstats$sim1[[1]][3,1,drop=TRUE]), network.dyadcount(nw))
  expect_equal(unname(sim$stats$nwstats$sim1[[1]][4,1,drop=TRUE]), 0)
})

test_that("tergm and tergmLite produce consistent durational statistics", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  nsims <- 3
  nsteps <- 5

  for (trim in list(FALSE, TRUE)) {
    if (trim == TRUE) {
      est2 <- trim_netest(est)
    } else {
      est2 <- est    
    }
    param <- param.net(inf.prob = 0, act.rate = 0)
    init <- init.net(i.num = 0)
    control <- control.net(type = NULL,
                           nsims = nsims,
                           nsteps = nsteps,
                           resimulate.network = TRUE,
                           tergmLite = FALSE,
                           nwstats.formula = ~edges + mean.age,
                           verbose = FALSE)

    set.seed(0)
    mod <- netsim(est2, param, init, control)
    stats <- get_nwstats(mod)

    control$tergmLite <- TRUE
    control$tergmLite.track.duration <- TRUE

    set.seed(0)
    mod2 <- netsim(est2, param, init, control)
    stats2 <- get_nwstats(mod2)

    expect_identical(stats, stats2)
    expect_true(all(stats$mean.age >= 0 & stats$mean.age <= nsteps + 1))

    control$tergmLite <- FALSE
    control$resimulate.network <- FALSE

    set.seed(0)
    mod3 <- netsim(est2, param, init, control)
    stats3 <- get_nwstats(mod3)
    expect_true(all(stats3$mean.age >= 0 & stats3$mean.age <= nsteps + 1))
  }
})

test_that("durational monitor with tergmLite", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

  nsims <- 2
  nsteps <- 4

  param <- param.net(inf.prob = 0.3, act.rate = 0.4)
  init <- init.net(i.num = 10)

  for (duration in list(1,10)) {
    est <- netest(nw, formation = ~edges + nodematch("race"),
                  target.stats = c(25, 10),
                  coef.diss = dissolution_coefs(~offset(edges), duration, 0),
                  verbose = FALSE)
    
    for (trim in list(FALSE, TRUE)) {
      if (trim == TRUE) {
        est2 <- trim_netest(est)
      } else {
        est2 <- est    
      }
      
      for (tergmLite.track.duration in list(FALSE, TRUE)) {
        control <- control.net(type = "SI",
                               nsims = nsims,
                               nsteps = nsteps,
                               resimulate.network = TRUE,
                               tergmLite = TRUE,
                               tergmLite.track.duration = tergmLite.track.duration,
                               nwstats.formula = ~edges + mean.age,
                               verbose = FALSE)
        
        mod <- netsim(est2, param, init, control)
        stats <- get_nwstats(mod)
        
        expect_true(all(stats$mean.age >= 0))
        print(stats$mean.age)
      }
    }
  }
})
