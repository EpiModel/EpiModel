context("New Network Models")

test_that("New network models vignette example", {
  skip_on_cran()

  ## New Aging Module
  aging <- function(dat, at) {

    age <- get_attr(dat, "age", override.null.error = TRUE)
    if (is.null(age)) {
      active <- get_attr(dat, "active")
      n <- sum(active == 1)
      age <- sample(18:49, n, replace = TRUE)
    } else {
      age <- get_attr(dat, "age") + 1 / 12
    }
    dat <- set_attr(dat, "age", age)

    return(dat)
  }


  ## Replacement Departure Module
  dfunc <- function(dat, at) {
    active <- get_attr(dat, "active")
    exitTime <- get_attr(dat, "exitTime")
    idsElig <- which(active == 1)
    nElig <- length(idsElig)

    nDepartures <- 0

    if (nElig > 0) {
      ages <- get_attr(dat, "age")[idsElig]
      life.expt <- get_param(dat, "life.expt")
      departure.rates <- pmin(1, 1 / (life.expt * 12 - ages * 12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        active[idsDepartures] <- 0
        exitTime[idsDepartures] <- at
        dat <- set_attr(dat, "active", active)
        dat <- set_attr(dat, "exitTime", exitTime)
      }
    }

    # Output
    dat <- set_epi(dat, "d.flow", at, nDepartures)
    return(dat)
  }


  ## Replacement Arrival Module
  afunc <- function(dat, at) {

    # Variables
    growth.rate <- get_param(dat, "growth.rate")
    exptPopSize <- get_epi(dat, "num", 1) * (1 + growth.rate * at)
    active <- get_attr(dat, "active")
    numNeeded <- exptPopSize - sum(active == 1)

    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }

    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 0, nArrivals)

    # Output
    dat <- set_epi(dat, "a.flow", at, nArrivals)

    return(dat)
  }


  ## Network Model
  nw <- network.initialize(50, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 15,
                coef.diss = dissolution_coefs(~offset(edges), 60, 0.000274),
                verbose = FALSE)


  ## EpiModel
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 5,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         tergmLite = FALSE, resimulate.network = TRUE, verbose = FALSE)
  mod1 <- netsim(est, param, init, control)
  capture_output(
    mod1
  )

  expect_is(mod1, "netsim")
  expect_output(print(mod1), "resim_nets.FUN")
  expect_output(print(mod1), "infection.FUN")
  expect_output(print(mod1), "departures.FUN")
  expect_output(print(mod1), "arrivals.FUN")
  expect_output(print(mod1), "aging.FUN")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = FALSE, resimulate.network = TRUE, verbose = FALSE)
  mod2 <- netsim(est, param, init, control)
  expect_is(mod2, "netsim")

  ### tergmLite replication
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         infection.FUN = infection.net,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         tergmLite = TRUE, verbose = FALSE,
                         resimulate.network = TRUE)
  mod3 <- netsim(est, param, init, control)
  expect_is(mod3, "netsim")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = TRUE, resimulate.network = TRUE, verbose = FALSE)
  mod4 <- netsim(est, param, init, control)
  expect_is(mod4, "netsim")

  ## "updated" infection module
  infect <- infection.net
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infect,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = TRUE, resimulate.network = TRUE, verbose = FALSE)
  mod5 <- netsim(est, param, init, control)
  expect_is(mod5, "netsim")

  expect_output(print(mod5), "resim_nets.FUN")
  expect_output(print(mod5), "infection.FUN")
  expect_output(print(mod5), "departures.FUN")
  expect_output(print(mod5), "arrivals.FUN")
  expect_output(print(mod5), "aging.FUN")

})

context("Network Model with Param Updater")

test_that("netsim with param updater", {
  skip_on_cran()
  # Create the list.param.updaters
  list.param.updaters <- list(
    # this is one updater
    list(
      at = 10,
      verbose = TRUE,
      param = list(
        inf.prob = 0.3,
        act.rate = 0.3
      )
    ),
    # this is another updater
    list(
      at = 20,
      verbose = TRUE,
      param = list(
        # inf.prob = function(x) plogis(qlogis(x) - log(10)),
        # act.rate = function(x) plogis(qlogis(x) - log(10))
        inf.prob = 0.01
      )
    )
  )

  # Create the list.control.updaters
  list.control.updaters <- list(
    # this is one updater
    list(
      at = 30,
      verbose = TRUE,
      control = list(
        resimulate.network = FALSE
      )
    )
  )

  # Do not forget to add it to `param`
  param <- param.net(
    inf.prob = 0.1,
    act.rate = 0.1,
    .param.updater.list = list.param.updaters
  )

  # Enable the module in `control`
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    .control.updater.list = list.control.updaters,
    resimulate.network = TRUE
  )

  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 25,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )

  init <- init.net(i.num = 10)

  expect_message(mod <- netsim(est, param, init, control))

  # `resimulate.network` is turned of at step 30. We check that the number of
  # observations in the "networkDynamic" object is < than 31 and not 50 (the
  # number of timestep in the simulation)
  n_obs <- length(
    get.network.attribute(mod$network[[1]][[1]], 'net.obs.period')$observations
  )
  expect_lt(n_obs, 31)
})

context("Network Model with Scenarios")

test_that("SIS with scenarios", {
  skip_on_cran()
  set.seed(10)

  nw <- network_initialize(n = 200)
  est <- netest(nw,
    formation = ~edges, target.stats = 60,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )

  param <- param.net(inf.prob = 0.9, rec.rate = 0.01, act.rate = 2)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 50, verbose = FALSE)
  init <- init.net(i.num = 10)

  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate,
    "base", 0, 0.9, 0.01,
    "multiple_changes", 0, 0.1, 0.04,
    "multiple_changes", 20, 0.9, 0.01,
    "multiple_changes", 40, 0.1, 0.1
  )

  scenarios.list <- create_scenario_list(scenarios.df)
  expect_length(scenarios.list, 2)

  sc.param <- use_scenario(param, scenarios.list[[1]])
  expect_silent(netsim(est, sc.param, init, control))
  sc.param <- use_scenario(param, scenarios.list[[2]])
  expect_message(netsim(est, sc.param, init, control))

  # .at not a integer
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate,
    "multiple_changes", "text", 0.1, 0.1
  )
  expect_error(create_scenario_list(scenarios.df))

  # inf_prob with an underscore
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf_prob, ~rec.rate,
    "multiple_changes", 0, 0.1, 0.1
  )
  expect_error(scenarios.list <- create_scenario_list(scenarios.df))

  # rec.rate2 not in param
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate2,
    "multiple_changes", 0, 0.1, 0.1
  )
  scenarios.list <- create_scenario_list(scenarios.df)
  expect_error(sc.param <- use_scenario(param, scenarios.list[[1]]))
})

context("Records: attr_history and Raw Objects")

test_that("Time varying elements", {
  skip_on_cran()
  test_logger <- function(dat, at) {
    nodes <- get_posit_ids(dat)

    some_nodes <- sample(nodes, 5)
    dat <- record_attr_history(
      dat, "attr_norm",
      rnorm(length(some_nodes)),
      posit_ids = some_nodes
    )

    some_nodes <- sample(nodes, 5)
    dat <- record_attr_history(
      dat, "attr_unif",
      runif(length(some_nodes)),
      posit_ids = some_nodes
    )

    some_nodes <- sample(nodes, 5)
    dat <- record_attr_history(
      dat, "attr_fix",
      at,
      posit_ids = some_nodes
    )

    # test when 0 nodes selected
    some_nodes <- integer(0)
    dat <- record_attr_history(
      dat, "attr_none",
      at,
      posit_ids = some_nodes
    )

    return(dat)
  }

  param <- param.net(
    inf.prob = 0.1,
    act.rate = 0.1
  )

  # Enable the module in `control`
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 20,
    verbose = FALSE,
    infection.FUN = infection.net,
    logger.FUN = test_logger
  )

  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 25,
    coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
    verbose = FALSE
  )

  init <- init.net(i.num = 10)

  mod <- netsim(est, param, init, control)
  attr_history <- get_attr_history(mod)
  expect_is(attr_history, "list")
  expect_is(attr_history[[1]], "data.frame")
  expect_equal(
    names(attr_history),
    c("attr_norm", "attr_unif", "attr_fix", "attr_none"))
})

context("Custom Trackers")

test_that("netsim, SI, custom trackers", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 25,
    coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
    verbose = FALSE
  )

  init <- init.net(i.num = 10)

  epi_s_num <- function(dat, at) {
    needed_attributes <- c("status")
    output <- with(get_attr_list(dat, needed_attributes), {
      sum(status == "s", na.rm = TRUE)
    })
    return(output)
  }

  epi_prop_infected <- function(dat, at) {
    needed_attributes <- c("status", "active")
    output <- with(get_attr_list(dat, needed_attributes), {
      pop <- active == 1
      cond <- status == "i"
      sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
    })
    return(output)
  }

  some.trackers <- list(
    prop_infected = epi_prop_infected,
    s_num = epi_s_num
  )

  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    .tracker.list = some.trackers
  )

  param <- param.net(
    inf.prob = 0.3,
    act.rate = 0.1
  )

  mod <- netsim(est, param, init, control)

  d <- as.data.frame(mod)

  expect_true(all(c("prop_infected", "s_num") %in% names(d)))
  expect_is(d[["prop_infected"]], "numeric")
  # the custom epi trackers are not run during intialization so the first value
  # is always NA
  expect_true(all(d$s_num[2:50] == d$s.num[2:50]))
})

context("Load Parameters from data.frame")

test_that("Load parameters from data.frame", {
  skip_on_cran()
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    "p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p3_1", "1", "numeric", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )

  expect_silent(param <- param.net_from_table(params.df))
  expect_s3_class(param, "param.net")
  expect_type(param, "list")

  expect_silent(param <- param.net(data.frame.parameters = params.df))
  expect_s3_class(param, "param.net")
  expect_type(param, "list")

  # convert back to a `long.param.df`
  param.df_back <- param.net_from_table(params.df) |> param.net_to_table()
  expect_true(all(params.df[c("param", "value", "type")] == param.df_back))

  # wrong column name
  params.df <- dplyr::tribble(~name, ~value, ~type, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))
  params.df <- dplyr::tribble(~param, ~val, ~type, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))
  params.df <- dplyr::tribble(~param, ~value, ~class, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))

  # wrong "type" value
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    "p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p3_1", "1", "factor", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )
  expect_error(param <- param.net_from_table(params.df))

  # wrong "param" format
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    ".p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p_3_1", "1", "numeric", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )
  expect_error(param <- param.net_from_table(params.df))
})

context("Random Parameter Generators")

test_that("Random parameters generators", {
  skip_on_cran()

  my_randoms <- list(
    act.rate = param_random(c(0.25, 0.5, 0.75)),
    tx.halt.part.prob = function() rbeta(1, 1, 2),
    hiv.test.rate = function() c(
      rnorm(1, 0.015, 0.01),
      rnorm(1, 0.010, 0.01),
      rnorm(1, 0.020, 0.01)
    )
  )

  expect_warning(param <- param.net(
      inf.prob = 0.3,
      act.rate = 0.3,
      random.params = my_randoms)
  )
  expect_message(generate_random_params(param, verbose = TRUE))
  expect_silent(generate_random_params(param, verbose = FALSE))

  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = list())
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = 4)
  expect_error(generate_random_params(param))

  param <- param.net(inf.prob = 0.3, random.params = list(1))
  expect_error(generate_random_params(param))


  generate_correlated_params <- function() {
    param.unique <- runif(1)
    param.set.1 <- param.unique + runif(2)
    param.set.2 <- param.unique * rnorm(3)

    return(list(param.unique, param.set.1, param.set.2))
  }

  # Data.frame set of random parameters :
  correlated_params <- t(replicate(10, unlist(generate_correlated_params())))
  correlated_params <- as.data.frame(correlated_params)
  colnames(correlated_params) <- c(
    "param.unique",
    "param.set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )

  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_silent(generate_random_params(param))

  # duplicated `act.rate` random definition
  colnames(correlated_params) <- c(
    "act.rate",
    "param.set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )
  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  expect_warning(
    param <- param.net(inf.prob = 0.3, act.rate = 0.1, random.params = randoms)
  )
  expect_warning(generate_random_params(param))

  # malformed name "param_set.1_1"
  colnames(correlated_params) <- c(
    "act.rate",
    "param_set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )
  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_error(generate_random_params(param))

  # param.random.set not a data.frame
  randoms <- c(my_randoms, list(param.random.set = list()))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_error(generate_random_params(param))
})
