test_that("mcmc controls default appropriately when NULL", {
  require(EpiModel)

  if (packageVersion("EpiModel") >= "2.1.0") {

    nw <- network_initialize(n = 50)
    nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

    ## tergm case
    est <- netest(nw, formation = ~edges + nodematch("race"),
                  target.stats = c(25, 10),
                  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                  verbose = FALSE)
    param <- param.net(inf.prob = 0.3, act.rate = 0.5)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI",
                           nsims = 1,
                           nsteps = 5,
                           verbose = FALSE,
                           tergmLite = TRUE,
                           resimulate.network = TRUE)
    set.seed(0)
    mod1 <- netsim(est, param, init, control)

    control <- control.net(type = "SI",
                           nsims = 1,
                           nsteps = 5,
                           verbose = FALSE,
                           tergmLite = TRUE,
                           resimulate.network = TRUE,
                           mcmc.control.tergm = NULL)
    set.seed(0)
    mod2 <- netsim(est, param, init, control)

    # should be equal except for the NULL control argument
    mod2$control$mcmc.control.tergm <- control.simulate.formula.tergm()

    expect_equal(mod1, mod2)

    ## ergm case
    est <- netest(nw, formation = ~edges + nodematch("race"),
                  target.stats = c(25, 10),
                  coef.diss = dissolution_coefs(~offset(edges), 1, 0),
                  verbose = FALSE)
    control <- control.net(type = "SI",
                           nsims = 1,
                           nsteps = 5,
                           verbose = FALSE,
                           tergmLite = TRUE,
                           resimulate.network = TRUE)
    set.seed(0)
    mod1 <- netsim(est, param, init, control)

    control <- control.net(type = "SI",
                           nsims = 1,
                           nsteps = 5,
                           verbose = FALSE,
                           tergmLite = TRUE,
                           resimulate.network = TRUE,
                           mcmc.control.ergm = NULL)
    set.seed(0)
    mod2 <- netsim(est, param, init, control)

    # should be equal except for the NULL control argument
    mod2$control$mcmc.control.ergm <- control.simulate.formula()

    expect_equal(mod1, mod2)

  }

})
