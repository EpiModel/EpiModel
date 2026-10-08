context("Tracked Attributes")

# -- Shared helpers ----------------------------------------------------------

simple_est <- function(n = 50) {
  nw <- network_initialize(n = n)
  netest(nw, formation = ~edges, target.stats = n / 2,
         coef.diss = dissolution_coefs(~offset(edges), 10, 0),
         verbose = FALSE)
}

# -- Tests -------------------------------------------------------------------

test_that("tracked.attributes records status changes in SI model", {
  skip_on_cran()
  est <- simple_est()
  param <- param.net(inf.prob = 0.3, act.rate = 2)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 20,
                         tracked.attributes = c("status"),
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)
  attr_hist <- get_attr_history(sim)

  # Both "active" and "status" should be tracked
  expect_true("active" %in% names(attr_hist))
  expect_true("status" %in% names(attr_hist))
  expect_is(attr_hist$status, "data.frame")

  # Initial snapshot: all 50 nodes recorded at first record
  init_records <- attr_hist$status[attr_hist$status$time == min(attr_hist$status$time), ]
  expect_equal(nrow(init_records), 50)

  # Subsequent steps should only have changes (fewer rows than 50)
  later_records <- attr_hist$status[attr_hist$status$time > min(attr_hist$status$time), ]
  if (nrow(later_records) > 0) {
    expect_true(all(later_records$values == "i"))  # SI: only s -> i
  }
})

test_that("tracked.attributes handles arrivals and departures", {
  skip_on_cran()

  aging <- function(dat, at) {
    age <- get_attr(dat, "age", override.null.error = TRUE)
    if (is.null(age)) {
      n <- sum(get_attr(dat, "active") == 1)
      age <- sample(18:49, n, replace = TRUE)
    } else {
      age <- get_attr(dat, "age") + 1 / 12
    }
    dat <- set_attr(dat, "age", age)
    return(dat)
  }

  dfunc <- function(dat, at) {
    active <- get_attr(dat, "active")
    exitTime <- get_attr(dat, "exitTime")
    idsElig <- which(active == 1)
    nDepartures <- 0
    if (length(idsElig) > 0) {
      ages <- get_attr(dat, "age")[idsElig]
      life.expt <- get_param(dat, "life.expt")
      departure.rates <- pmin(1, 1 / (life.expt * 12 - ages * 12))
      vecDep <- which(rbinom(length(idsElig), 1, departure.rates) == 1)
      idsDep <- idsElig[vecDep]
      nDepartures <- length(idsDep)
      if (nDepartures > 0) {
        active[idsDep] <- 0
        exitTime[idsDep] <- at
        dat <- set_attr(dat, "active", active)
        dat <- set_attr(dat, "exitTime", exitTime)
      }
    }
    dat <- set_epi(dat, "d.flow", nDepartures)
    return(dat)
  }

  afunc <- function(dat, at) {
    growth.rate <- get_param(dat, "growth.rate")
    exptPopSize <- get_epi(dat, "num", 1) * (1 + growth.rate * at)
    numNeeded <- exptPopSize - sum(get_attr(dat, "active") == 1)
    nArrivals <- if (numNeeded > 0) rpois(1, numNeeded) else 0
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 0, nArrivals)
    dat <- set_epi(dat, "a.flow", nArrivals)
    return(dat)
  }

  nw <- network.initialize(50, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 15,
                coef.diss = dissolution_coefs(~offset(edges), 60, 0.000274),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 20,
                         departures.FUN = dfunc, arrivals.FUN = afunc,
                         aging.FUN = aging, infection.FUN = infection.net,
                         tergmLite = FALSE, resimulate.network = TRUE,
                         tracked.attributes = c("status"),
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)
  attr_hist <- get_attr_history(sim)
  d_active <- attr_hist$active

  # Initial snapshot: all nodes active = 1
  init_active <- d_active[d_active$time == min(d_active$time), ]
  expect_true(all(init_active$values == 1L))

  # Departures recorded as active = 0
  departures <- d_active[d_active$time > min(d_active$time) & d_active$values == 0L, ]
  arrivals <- d_active[d_active$time > min(d_active$time) & d_active$values == 1L, ]

  # Arrivals should have UIDs > initial population
  if (nrow(arrivals) > 0) {
    expect_true(all(arrivals$uids > 50))
  }

  # Arrivals should also appear in status history
  if (nrow(arrivals) > 0) {
    arrival_uids <- unique(arrivals$uids)
    status_uids <- unique(attr_hist$status$uids)
    expect_true(all(arrival_uids %in% status_uids))
  }
})

test_that("get_attr_at reconstructs correct state", {
  skip_on_cran()
  est <- simple_est()
  param <- param.net(inf.prob = 0.3, act.rate = 2)
  init <- init.net(i.num = 10)
  nsteps <- 20
  control <- control.net(type = "SI", nsims = 1, nsteps = nsteps,
                         tracked.attributes = c("status"),
                         save.other = c("attr"),
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)

  # Reconstruct at final step
  attrs_end <- get_attr_at(sim, at = nsteps)

  # Compare against saved attributes
  final_attr <- sim$attr[[1]]
  active_ids <- which(final_attr$active == 1)
  final_uids <- final_attr$unique_id[active_ids]
  final_status <- final_attr$status[active_ids]

  expect_equal(sort(attrs_end$unique_id), sort(final_uids))
  reconstructed <- attrs_end$status[match(final_uids, attrs_end$unique_id)]
  expect_equal(reconstructed, final_status)
})

test_that("get_attr_at validates inputs", {
  skip_on_cran()
  est <- simple_est()
  param <- param.net(inf.prob = 0.3, act.rate = 2)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 10,
                         tracked.attributes = c("status"),
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_error(get_attr_at(sim, at = 0))
  expect_error(get_attr_at(sim, at = 11))
  expect_error(get_attr_at(sim, at = 5, sim_num = 2))
  expect_error(get_attr_at("not_a_sim", at = 5))

  # Without active tracking
  control2 <- control.net(type = "SI", nsims = 1, nsteps = 10,
                          verbose = FALSE)
  sim2 <- netsim(est, param, init, control2)
  expect_error(get_attr_at(sim2, at = 5), "tracked")
})

test_that("empty tracked.attributes does not create attr.history", {
  skip_on_cran()
  est <- simple_est()
  param <- param.net(inf.prob = 0.3, act.rate = 2)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 10,
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)
  attr_hist <- get_attr_history(sim)
  expect_equal(length(attr_hist), 0)
})

test_that("tracked.attributes and manual recording coexist", {
  skip_on_cran()

  manual_logger <- function(dat, at) {
    dat <- record_attr_history(dat, "custom_metric", at * 2,
                               unique_ids = get_unique_ids(dat)[1])
    return(dat)
  }

  est <- simple_est()
  param <- param.net(inf.prob = 0.3, act.rate = 2)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         infection.FUN = infection.net,
                         logger.FUN = manual_logger,
                         tracked.attributes = c("status"),
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)
  attr_hist <- get_attr_history(sim)

  # All three should be present
  expect_true("active" %in% names(attr_hist))
  expect_true("status" %in% names(attr_hist))
  expect_true("custom_metric" %in% names(attr_hist))
})
