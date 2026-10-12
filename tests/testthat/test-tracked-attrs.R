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

test_that("get_attr_at computes age only with a time.unit parameter", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "age", runif(50, 20, 40))
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  aging <- function(dat, at) {
    dat <- set_attr(dat, "age", get_attr(dat, "age") + 1 / 12)
    return(dat)
  }
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 20,
                         aging.FUN = aging, infection.FUN = infection.net,
                         tracked.attributes.once = "age",
                         save.other = "attr", verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 2)
  expect_null(param$time.unit)
  sim <- netsim(est, param, init, control)
  expect_error(get_attr_at(sim, at = 20, compute_age = TRUE), "time.unit")
  # Without compute_age, age is the last recorded value
  attrs <- get_attr_at(sim, at = 20)
  first_age <- get_attr_history(sim)$age
  expect_equal(attrs$age, first_age$values[match(attrs$unique_id, first_age$uids)])

  # Monthly time steps
  param <- param.net(inf.prob = 0.3, act.rate = 2, time.unit = 364 / 12)
  sim <- netsim(est, param, init, control)
  attrs <- get_attr_at(sim, at = 20, compute_age = TRUE)
  final <- sim$attr[[1]]
  expect_equal(attrs$age, final$age[match(attrs$unique_id, final$unique_id)])
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

# -- Reconstruction against the true attribute list --------------------------

# An open-population SI model that saves the attribute list at the end of each
# step with `record_raw_object()`, as the truth to compare the history against
open_pop_fixture <- function(n = 100, nsteps = 30,
                             tracked = c("status", "infTime"),
                             once = "entrTime", ...) {
  set.seed(42)
  nw <- network_initialize(n)
  est <- netest(nw, ~edges, n / 2, dissolution_coefs(~offset(edges), 10, 0.02),
                verbose = FALSE)
  prev_log <- function(dat, at) {
    dat <- prevalence.net(dat, at)
    record_raw_object(
      dat, "attrs",
      get_attr_list(dat, c("unique_id", "active", "status", "infTime"))
    )
  }
  param <- param.net(inf.prob = 0.3, act.rate = 1, a.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = nsteps,
                         infection.FUN = infection.net,
                         departures.FUN = departures.net,
                         arrivals.FUN = arrivals.net,
                         prevalence.FUN = prev_log,
                         resimulate.network = TRUE, tergmLite = TRUE,
                         save.run = TRUE,
                         tracked.attributes = tracked,
                         tracked.attributes.once = once, verbose = FALSE, ...)
  list(est = est, param = param, init = init, control = control,
       sim = netsim(est, param, init, control))
}

truth_at <- function(sim, at) {
  recs <- sim$raw.records$sim1$as_list()
  a <- recs[[which(vapply(recs, `[[`, 0, "time") == at)]]$object
  keep <- a$active == 1
  data.frame(unique_id = a$unique_id[keep], status = a$status[keep],
             infTime = a$infTime[keep])
}

expect_matches_truth <- function(sim, at) {
  got <- get_attr_at(sim, at)
  want <- truth_at(sim, at)
  expect_setequal(got$unique_id, want$unique_id)
  m <- match(want$unique_id, got$unique_id)
  expect_identical(got$status[m], want$status)
  expect_equal(got$infTime[m], want$infTime)
}

test_that("get_attr_at matches the attribute list at every step", {
  skip_on_cran()
  f <- open_pop_fixture()
  expect_true(any(f$sim$epi$a.flow[, 1] > 0, na.rm = TRUE))
  expect_true(any(f$sim$epi$ds.flow[, 1] + f$sim$epi$di.flow[, 1] > 0,
                  na.rm = TRUE))
  for (at in 2:30) {
    expect_matches_truth(f$sim, at)
  }

  # the initial snapshot is recorded at start - 1, arrivals at their entrTime
  h <- get_attr_history(f$sim)
  expect_equal(min(h$active$time), 1)
  expect_equal(sum(h$active$time == 1), 100)
  arr <- h$active[h$active$time > 1 & h$active$values == 1, ]
  ent <- h$entrTime[match(arr$uids, h$entrTime$uids), ]
  expect_equal(ent$values, arr$time)
})

test_that("tracking starts again from a full snapshot on restart", {
  skip_on_cran()
  f <- open_pop_fixture(nsteps = 10)
  ctrl2 <- f$control
  ctrl2$start <- 11
  ctrl2$nsteps <- 20
  sim2 <- netsim(f$sim, f$param, f$init, ctrl2)
  h <- get_attr_history(sim2)
  expect_equal(min(h$active$time), 10)
  for (at in c(11, 15, 20)) {
    expect_matches_truth(sim2, at)
  }
  # the restart does not add records to the original object
  expect_equal(max(get_attr_history(f$sim)$active$time), 10)

  rp <- make_restart_point(f$sim, time_attrs = "infTime")
  expect_false(rp$run[[1]]$tracking_attrs)
  expect_type(rp$run[[1]]$attr$unique_id, "integer")
  ctrl3 <- f$control
  ctrl3$start <- 2
  ctrl3$nsteps <- 11
  sim3 <- netsim(rp, f$param, f$init, ctrl3)
  expect_matches_truth(sim3, 11)
})

test_that("get_nodes_spell matches unique IDs by value", {
  # 1e5 prints as "1e+05", so a lookup by name never found its departure
  d <- data.frame(time = c(1, 1, 1, 5, 5, 5),
                  uids = c(99999, 1e5, 100001, 99999, 1e5, 100001),
                  values = c(1, 1, 1, 0, 0, 0))
  sp <- EpiModel:::get_nodes_spell(d)
  expect_equal(sp$uid, c(99999, 1e5, 100001))
  expect_equal(sp$terminus, c(5, 5, 5))

  # a node that is active again gets a second spell, and a repeated departure
  # does not end a spell twice
  d <- data.frame(time = c(1, 1, 3, 3, 5, 6, 7),
                  uids = c(1, 2, 1, 2, 1, 2, 2),
                  values = c(1, 1, 0, 0, 1, 0, 0))
  sp <- EpiModel:::get_nodes_spell(d)
  expect_equal(sp$uid, c(1, 1, 2))
  expect_equal(sp$onset, c(1, 5, 1))
  expect_equal(sp$terminus, c(3, Inf, 3))
  expect_equal(EpiModel:::get_spell_terminus(sp, c(1, 1, 2), c(2, 6, 1)),
               c(3, Inf, 3))
})

test_that("tracked.attributes.once alone tracks active and arrivals", {
  skip_on_cran()
  f <- open_pop_fixture(tracked = character(0), once = "entrTime")
  h <- get_attr_history(f$sim)
  expect_named(h, c("active", "entrTime"))
  expect_setequal(get_attr_at(f$sim, 30)$unique_id,
                  truth_at(f$sim, 30)$unique_id)
})

test_that("tracked attributes are checked", {
  expect_error(control.net(type = "SI", nsteps = 5, tracked.attributes = 1),
               "must be a character vector")
  expect_error(control.net(type = "SI", nsteps = 5,
                           tracked.attributes.once = c("race", NA)),
               "must be a character vector")
  expect_error(control.net(type = "SI", nsteps = 5,
                           tracked.attributes = c("status", "unique_id")),
               "cannot include \"unique_id\"")
  skip_on_cran()
  expect_error(open_pop_fixture(nsteps = 3, tracked = "statsu"),
               "do not exist at the end of time step 2: `statsu`")
})

test_that("attributes created by a module on its first call are tracked", {
  skip_on_cran()
  set.seed(1)
  est <- simple_est(30)
  mk <- function(dat, at) {
    n <- length(get_attr(dat, "active"))
    if (at == 2) {
      dat <- set_attr(dat, "late", rep(1, n))
      dat <- set_attr(dat, "late_once", seq_len(n))
    }
    if (at == 4) {
      late <- get_attr(dat, "late")
      late[1:3] <- 2
      dat <- set_attr(dat, "late", late)
    }
    dat
  }
  control <- control.net(type = NULL, nsims = 1, nsteps = 5, late.FUN = mk,
                         infection.FUN = infection.net,
                         tracked.attributes = c("status", "late"),
                         tracked.attributes.once = "late_once",
                         verbose = FALSE)
  sim <- netsim(est, param.net(inf.prob = 0.1, act.rate = 1),
                init.net(i.num = 1), control)
  h <- get_attr_history(sim)
  expect_named(h, c("active", "status", "late", "late_once"))
  expect_equal(unique(h$late$time), c(2, 4))
  expect_equal(sum(h$late$time == 2), 30)
  expect_equal(h$late$uids[h$late$time == 4], 1:3)
  expect_true(all(is.na(get_attr_at(sim, 1)$late)))
  expect_equal(get_attr_at(sim, 5)$late, rep(2:1, c(3, 27)))
  expect_equal(get_attr_at(sim, 5)$late_once, 1:30)
})

test_that("get_attr_history keeps the order of the first records", {
  skip_on_cran()
  logger <- function(dat, at) {
    dat <- record_attr_history(dat, "zz", 1, unique_ids = 1)
    dat <- record_attr_history(dat, "aa", 2, unique_ids = 1)
    return(dat)
  }
  control <- control.net(type = NULL, nsims = 1, nsteps = 3,
                         infection.FUN = infection.net, logger.FUN = logger,
                         tracked.attributes = "status",
                         tracked.attributes.once = "entrTime",
                         verbose = FALSE)
  sim <- netsim(simple_est(), param.net(inf.prob = 0.3, act.rate = 2),
                init.net(i.num = 10), control)
  expect_named(get_attr_history(sim),
               c("active", "status", "entrTime", "zz", "aa"))
})

test_that("truncate_sim keeps get_attr_at aligned with the new time steps", {
  skip_on_cran()
  f <- open_pop_fixture(nsteps = 20, cumulative.edgelist = TRUE,
                        truncate.el.cuml = Inf,
                        save.cumulative.edgelist = TRUE)
  tr <- truncate_sim(f$sim, at = 11)
  expect_identical(get_attr_at(tr, 1), get_attr_at(f$sim, 11))
  expect_identical(get_attr_at(tr, 5), get_attr_at(f$sim, 15))
  expect_identical(get_attr_at(tr, 10), get_attr_at(f$sim, 20))
  el <- f$sim$cumulative.edgelist$sim1
  el_tr <- tr$cumulative.edgelist$sim1
  expect_equal(el_tr$start, el$start - 10)
  expect_equal(el_tr$stop, el$stop - 10)
  # the original object is left as is
  expect_equal(min(get_attr_history(f$sim)$active$time), 1)
  expect_equal(min(get_attr_history(tr)$active$time), -9)
  expect_identical(truncate_sim(f$sim, at = 11, reset.time = FALSE)$attr.history,
                   f$sim$attr.history)
})

test_that("get_attr_at errors clearly without an attribute history", {
  skip_on_cran()
  f <- open_pop_fixture(nsteps = 5)
  sim <- f$sim
  sim$attr.history <- NULL
  expect_error(get_attr_at(sim, 5), "no attribute history")
  expect_error(get_attr_at(f$sim, c(2, 3)), "single at")
})

test_that("make_networkDynamic refuses a truncated cumulative edgelist", {
  skip_on_cran()
  set.seed(1)
  est <- simple_est()
  control <- control.net(type = "SI", nsims = 1, nsteps = 10,
                         cumulative.edgelist = TRUE,
                         save.cumulative.edgelist = TRUE,
                         tracked.attributes = "status", verbose = FALSE)
  sim <- netsim(est, param.net(inf.prob = 0.1, act.rate = 1),
                init.net(i.num = 5), control)
  expect_error(make_networkDynamic(sim), "truncate.el.cuml = Inf")
})
