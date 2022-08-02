context("Standard network models")

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(nw, formation = ~edges + nodematch("race"),
              target.stats = c(25, 10),
              coef.diss = dissolution_coefs(~offset(edges), 10, 0),
              verbose = FALSE)

for (trim in c(FALSE, TRUE)) {
  if (trim == TRUE) {
    est2 <- trim_netest(est)
  } else {
    est2 <- est
  }  
  
  # Edges + nodematch, one-mode, closed
  
  test_that("netsim for edges only, SI, one-mode, closed, 1 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
    plot(mod)
    plot(mod, type = "formation")
    plot(mod, type = "duration")
    plot(mod, type = "dissolution")    
    plot(mod, type = "network")
    test_net(mod)
  })
  
  test_that("netsim for edges only, SI, one-mode, closed, 2 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
    plot(mod)
    plot(mod, type = "formation")
    plot(mod, type = "network")
    test_net(mod)
  })
  
  test_that("netsim for edges only, SIS, one-mode, closed, 1 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SIS", nsims = 1, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
    plot(mod)
    plot(mod, type = "formation")
    plot(mod, type = "network")
    test_net(mod)
  })
  
  test_that("netsim for edges only, SIS, one-mode, closed, 2 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SIS", nsims = 2, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
  })
  
  test_that("netsim for edges only, SIR, one-mode, closed, 1 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
    init <- init.net(i.num = 10, r.num = 0)
    control <- control.net(type = "SIR", nsims = 1, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
    plot(mod)
    plot(mod, type = "formation")
    plot(mod, type = "network")
    test_net(mod)
  })
  
  test_that("netsim for edges only, SIR, one-mode, closed, 2 sim", {
    param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
    init <- init.net(i.num = 10, r.num = 0)
    control <- control.net(type = "SIR", nsims = 2, nsteps = 5, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
    plot(mod)
    plot(mod, type = "formation")
    plot(mod, type = "network")
    test_net(mod)
  })
}

test_that("netsim implicit save.network option", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  mod1 <- netsim(est1, param, init, control)
  expect_s3_class(get_network(mod1), "networkDynamic")

  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE,
                         save.network = FALSE)
  mod2 <- netsim(est1, param, init, control)
  expect_error(get_network(mod2))
})

test_that("netsim duration 1", {
  estd1 <- netest(nw, formation = ~edges + nodematch("race"),
                  target.stats = c(25, 10),
                  coef.diss = dissolution_coefs(~offset(edges), 1, 0),
                  verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(r.num = 0, status.vector = rep("s", 50))
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, 
                         resimulate.network = TRUE, verbose = FALSE,
                         nwupdate.FUN = NULL)
  set.seed(0)
  mod <- netsim(estd1, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
  
  # compare to manually produced networkDynamic
  set.seed(0)
  sim <- simulate(estd1$formation,
                  coef = estd1$coef.form.crude,
                  basis = estd1$newnetwork,
                  control = control.simulate.formula(MCMC.burnin = 2e5),
                  dynamic = FALSE)
  for(i in 1:5) {
    sim <- simulate(estd1$formation,
                    basis = sim,
                    time.slices = 1,
                    time.start = i - 1,
                    time.offset = 1,
                    coef = c(estd1$coef.form),
                    dynamic = TRUE)
  }
  expect_identical(sim$mel, mod$network$sim1[[1]]$mel)  
})

test_that("non-nested EDA works in netsim", {
  nw <- network.initialize(10, directed = FALSE)
  nw %v% "race" <- rep(1:2, length.out = 10)
  nw %v% "age" <- rep(1:5, length.out = 10)
  dc <- dissolution_coefs(~offset(edges) + offset(nodematch("age")), c(3, 7))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(10, 5),
                coef.diss = dc, nested.edapprox = FALSE)
  dxs <- netdx(est, nsteps = 2, nsims = 2, dynamic = FALSE)
  dxd <- netdx(est, nsteps = 2, nsims = 2, dynamic = TRUE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  
  dc <- dissolution_coefs(~offset(edges) + offset(nodematch("age")), c(1, 1))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(10, 5),
                coef.diss = dc, nested.edapprox = FALSE)
  dxs <- netdx(est, nsteps = 2, nsims = 2, dynamic = FALSE)
  sim <- netsim(est, param, init, control)
})

test_that("netsim diss.stats", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  mod <- netsim(est, param, init, control)

  print(mod)
  expect_output(print(mod), "Duration Statistics.*Sim Mean")
  expect_output(print(mod), "Dissolution Statistics.*Sim Mean")
  expect_error(expect_output(print(mod), "Not available when:"))
  expect_equal(length(mod[["diss.stats"]]), 2)
  
  mod2 <- merge(mod, mod)
  expect_output(print(mod2), "Duration Statistics.*Sim Mean")
  expect_output(print(mod2), "Dissolution Statistics.*Sim Mean")
  expect_error(expect_output(print(mod2), "Not available when:"))
  expect_equal(length(mod2[["diss.stats"]]), 4)
  
  mod3 <- merge(mod, mod, keep.diss.stats = FALSE)
  expect_output(print(mod3), "Duration and Dissolution Statistics")
  expect_output(print(mod3), "Not available when:")
  expect_true(is.null(mod3[["diss.stats"]]))
})

test_that("save.other sim naming", {
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges,
                target.stats = c(25),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE, 
                         save.other = c("nw"), resimulate.network = TRUE)
  mod <- netsim(est, param, init, control)
  expect_equal(names(mod[["nw"]]), paste0("sim", 1:2))
  mod2 <- merge(mod, mod)
  expect_equal(names(mod2[["nw"]]), paste0("sim", 1:4))
  mod3 <- get_sims(mod2, c(1,3,4))
  expect_equal(names(mod3[["nw"]]), paste0("sim", 1:3))
  mod4 <- merge(mod, mod, keep.other = FALSE)
  expect_equal(names(mod4[["nw"]]), NULL)
})
