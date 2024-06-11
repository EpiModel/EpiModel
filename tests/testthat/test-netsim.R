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
  sim <- networkDynamic::as.networkDynamic(sim)
  sim <- networkDynamic::activate.vertices(sim, onset = 0, terminus = Inf)
  sim <- networkDynamic::activate.edges(sim, onset = 0, terminus = Inf)

  sim %n% "net.obs.period" <- list(observations = list(c(0,1)),
                                   mode = "discrete",
                                   time.increment =  1,
                                   time.unit = "step")

  for(i in 1:5) {
    sim <- suppressWarnings(simulate(estd1$formation,
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i,
                                     time.offset = 0,
                                     control = list(MCMC.prop.args = list(discordance_fraction = 0)),
                                     coef = c(estd1$coef.form),
                                     dynamic = TRUE))
  }
  expect_identical(as.data.frame(sim), as.data.frame(mod$network$sim1[[1]]))
})

test_that("non-nested EDA works in netsim", {
  nw <- network.initialize(10, directed = FALSE)
  nw %v% "race" <- rep(1:2, length.out = 10)
  nw %v% "age" <- rep(1:5, length.out = 10)
  dc <- dissolution_coefs(~offset(edges) + offset(nodematch("age")), c(3, 7))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(10, 5),
                coef.diss = dc, nested.edapprox = FALSE)
  dxs <- netdx(est, nsteps = 2, nsims = 2, dynamic = FALSE, verbose = FALSE)
  dxd <- netdx(est, nsteps = 2, nsims = 2, dynamic = TRUE, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  sim <- netsim(est, param, init, control)

  dc <- dissolution_coefs(~offset(edges) + offset(nodematch("age")), c(1, 1))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(10, 5),
                coef.diss = dc, nested.edapprox = FALSE)
  dxs <- netdx(est, nsteps = 2, nsims = 2, dynamic = FALSE, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")
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

  capture_output(
    print(mod)
  )
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

test_that("name_saveout_elts unit", {
  simnames <- paste0("sim", 1:4)
  elt_name <- "this_elt"

  elt <- rep(list(sample(10)), 4)
  named_elt <- expect_silent(name_saveout_elts(elt, elt_name, simnames))
  expect_equal(names(named_elt), simnames)

  # wrong size produces a warning
  elt <- rep(list(sample(10)), 2)
  named_elt <- expect_warning(name_saveout_elts(elt, elt_name, simnames))
  expect_null(names(named_elt))
  expect_equal(elt, named_elt)

  # empty element returns silently
  elt <- NULL
  named_elt <- expect_silent(name_saveout_elts(elt, elt_name, simnames))
  expect_null(names(named_elt))
  expect_equal(elt, named_elt)
})

test_that("edges correction behaves as expected", {
  nw <- network_initialize(n = 500)

  est <- netest(nw, formation = ~edges + degree(1),
                target.stats = c(250, 100),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param_g1 <- param.net(inf.prob = 0.4, act.rate = 1,
                        a.rate = 0.03, ds.rate = 0.03, di.rate = 0.03)

  param_g2 <- param.net(inf.prob = 0.4, act.rate = 1, inf.prob.g2 = 0.4,
                        a.rate = 0.03, ds.rate = 0.03, di.rate = 0.03,
                        a.rate.g2 = 0.03, ds.rate.g2 = 0.03, di.rate.g2 = 0.03)

  nsims <- 2

  init_g1 <- init.net(i.num = 100)
  init_g2 <- init.net(i.num = 100, i.num.g2 = 100)
  for (tergmLite in list(FALSE, TRUE)) {
    for (resimulate.network in unique(c(tergmLite, TRUE))) {
      for (ngroups in list(1,2)) {
        if (ngroups == 1) {
          est$newnetwork <- delete.vertex.attribute(est$newnetwork, "group")
          param <- param_g1
          init <- init_g1
        } else {
          est$newnetwork %v% "group" <- rep(c(1,2), length.out = network.size(est$newnetwork))
          param <- param_g2
          init <- init_g2
        }

        nsteps <- 5
        control <- control.net(type = "SI", nsteps = nsteps, nsims = nsims, ncores = 1,
                               resimulate.network = resimulate.network,
                               tergmLite = tergmLite,
                               verbose = FALSE,
                               save.network = TRUE,
                               save.run = TRUE,
                               save.other = c("temp"))
        sim <- netsim(est, param, init, control)

        for (simno in seq_len(nsims)) {
          if (ngroups == 1) {
            expect_equal(est$coef.form[1] + log(network.size(nw)),
                         sim$coef.form[[simno]][[1]][1] + log(sim$run[[simno]]$num),
                         tolerance = 1e-6)
          } else {
            n1.old <- sum(est$newnetwork %v% "group" == 1)
            n2.old <- sum(est$newnetwork %v% "group" == 2)
            n1.new <- sim$run[[simno]]$num
            n2.new <- sim$run[[simno]]$num.g2

            expect_equal(est$coef.form[1] + log(2*n1.old*n2.old/(n1.old+n2.old)),
                         sim$coef.form[[simno]][[1]][1] + log(2*n1.new*n2.new/(n1.new+n2.new)),
                         tolerance = 1e-6)
          }
        }

        if (tergmLite == FALSE) {
          control$start <- nsteps + 1
          nsteps <- 9
          control$nsteps <- nsteps
          sim2 <- netsim(sim, param, init, control)

          for (simno in seq_len(nsims)) {
            if (ngroups == 1) {
              expect_equal(est$coef.form[1] + log(network.size(nw)),
                          sim$coef.form[[simno]][[1]][1] + log(sim$run[[simno]]$num),
                          tolerance = 1e-6)
            } else {
              n1.old <- sum(est$newnetwork %v% "group" == 1)
              n2.old <- sum(est$newnetwork %v% "group" == 2)
              n1.new <- sim$run[[simno]]$num
              n2.new <- sim$run[[simno]]$num.g2

              expect_equal(est$coef.form[1] + log(2*n1.old*n2.old/(n1.old+n2.old)),
                          sim$coef.form[[simno]][[1]][1] + log(2*n1.new*n2.new/(n1.new+n2.new)),
                          tolerance = 1e-6)
              }
          }
        }
      }
    }
  }
})

test_that("networkDynamics produced by netsim match those produced by simulate when they should", {
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges,
                target.stats = c(25),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0, act.rate = 0, rec.rate = 0)
  init <- init.net(i.num = 0, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, verbose = FALSE,
                         save.network = TRUE, resimulate.network = TRUE,
                         save.diss.stats = FALSE,
                         save.run = TRUE,
                         save.other = c("temp"),
                         save.transmat = FALSE)
  set.seed(0)
  mod <- netsim(est, param, init, control)
  control$start <- 6
  control$nsteps <- 11
  mod2 <- netsim(mod, param, init, control)

  # compare to manually produced networkDynamic
  set.seed(0)
  sim <- simulate(est$formation,
                  coef = est$coef.form.crude,
                  basis = est$newnetwork,
                  control = control.simulate.formula(MCMC.burnin = 2e5),
                  dynamic = FALSE)
  sim <- networkDynamic::as.networkDynamic(sim)
  sim <- networkDynamic::activate.vertices(sim, onset = 0, terminus = Inf)
  sim <- networkDynamic::activate.edges(sim, onset = 0, terminus = Inf)

  sim %n% "net.obs.period" <- list(observations = list(c(0,1)),
                                   mode = "discrete",
                                   time.increment =  1,
                                   time.unit = "step")

  for(i in 1:5) {
    sim <- suppressWarnings(simulate(~Form(est$formation) + Persist(est$coef.diss$dissolution),
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i,
                                     time.offset = 0,
                                     coef = c(est$coef.form, est$coef.diss$coef.crude),
                                     dynamic = TRUE))
  }
  expect_identical(as.data.frame(sim), as.data.frame(mod$network$sim1[[1]]))

  for(i in 6:11) {
    sim <- suppressWarnings(simulate(~Form(est$formation) + Persist(est$coef.diss$dissolution),
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i,
                                     time.offset = 0,
                                     coef = c(est$coef.form, est$coef.diss$coef.crude),
                                     dynamic = TRUE))
  }
  expect_identical(as.data.frame(sim), as.data.frame(mod2$network$sim1[[1]]))
})

test_that("networkLites produced by netsim match those produced by simulate when they should", {
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges,
                target.stats = c(25),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0, act.rate = 0, rec.rate = 0)
  init <- init.net(i.num = 0, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, verbose = FALSE,
                         save.network = TRUE, resimulate.network = TRUE,
                         tergmLite = TRUE,
                         save.run = TRUE,
                         save.other = c("temp", "el", "net_attr"),
                         tergmLite.track.duration = TRUE, save.transmat = FALSE)
  est <- trim_netest(est)
  set.seed(0)
  mod <- netsim(est, param, init, control)
  control$start <- 6
  control$nsteps <- 11
  mod2 <- netsim(mod, param, init, control)

  # compare to manually produced networkDynamic
  set.seed(0)
  sim <- simulate(est$formation,
                  coef = est$coef.form.crude,
                  basis = est$newnetwork,
                  control = control.simulate.formula(MCMC.burnin = 2e5),
                  dynamic = FALSE)

  sim %n% "time" <- 0L
  sim %n% "lasttoggle" <- cbind(as.edgelist(sim), 0L)

  for(i in 1:5) {
    sim <- suppressWarnings(simulate(~Form(est$formation) + Persist(est$coef.diss$dissolution),
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i - 1,
                                     time.offset = 1,
                                     output = "final",
                                     coef = c(est$coef.form, est$coef.diss$coef.crude),
                                     dynamic = TRUE))
  }
  expect_equal(as.edgelist(sim), as.edgelist(mod$network$sim1[[1]]), check.attributes = FALSE)
  expect_equal(sim %n% "lasttoggle", mod$network$sim1[[1]] %n% "lasttoggle")

  for(i in 6:11) {
    sim <- suppressWarnings(simulate(~Form(est$formation) + Persist(est$coef.diss$dissolution),
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i - 1,
                                     time.offset = 1,
                                     output = "final",
                                     coef = c(est$coef.form, est$coef.diss$coef.crude),
                                     dynamic = TRUE))
  }
  expect_equal(as.edgelist(sim), as.edgelist(mod2$network$sim1[[1]]), check.attributes = FALSE)
  expect_equal(sim %n% "lasttoggle", mod2$network$sim1[[1]] %n% "lasttoggle")
})
