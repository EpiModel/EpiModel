context("Observed Network Layers")

# Shared fixtures ------------------------------------------------------------

# A static census: a fixed observed network with a nodal attribute
set.seed(41)
N <- 60
nw_static <- network_initialize(n = N)
nw_static <- set_vertex_attribute(nw_static, "sex", rep(c("F", "M"), N / 2))
pairs <- t(combn(N, 2))
pairs <- pairs[sample.int(nrow(pairs), 90), ]
nw_static <- add.edges(nw_static, tail = pairs[, 1], head = pairs[, 2])

# A dynamic census: edge spells over an observation window of 0 to 30. Each
# observed dyad is active for one interval; some start before the window
# (censored onsets) and some are still active at its end.
make_dynamic_census <- function(n = 60, n.edges = 150, end = 30) {
  pairs <- t(combn(n, 2))
  pairs <- pairs[sample.int(nrow(pairs), n.edges), ]
  onset <- sample(c(-Inf, 0:(end - 5)), n.edges, replace = TRUE)
  terminus <- pmax(ifelse(is.finite(onset), onset, 0), 0) +
    sample(2:12, n.edges, replace = TRUE)
  terminus[terminus >= end] <- Inf
  base <- network::network.initialize(n, directed = FALSE)
  base <- set_vertex_attribute(base, "sex", rep(c("F", "M"), n / 2))
  nd <- networkDynamic::networkDynamic(
    base.net = base,
    edge.spells = data.frame(onset = onset, terminus = terminus,
                             tail = pairs[, 1], head = pairs[, 2])
  )
  nd %n% "net.obs.period" <- list(observations = list(c(0, end)),
                                  mode = "discrete", time.increment = 1,
                                  time.unit = "step")
  nd
}
nw_dyn <- suppressMessages(make_dynamic_census())

active_at <- function(nd, at) {
  el <- networkDynamic::get.dyads.active(nd, at = at)
  el <- cbind(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]))
  el[order(el[, 1], el[, 2]), , drop = FALSE]
}

init <- init.net(i.num = 5)


# Constructor ------------------------------------------------------------------

test_that("netcensus wraps a static network", {
  obs <- netcensus(nw_static)
  expect_s3_class(obs, "netcensus")
  expect_false(obs$dynamic)
  expect_null(obs$window)
  expect_false(obs$edapprox)
  expect_null(obs$target.stats)
  expect_equal(obs$summary$edges, 90)
  expect_equal(obs$summary$n, N)
  expect_equal(obs$summary$mean.degree, 2 * 90 / N)
  expect_equal(obs$summary$attributes, "sex")
  out <- capture.output(print(obs))
  expect_true(any(grepl("static \\(network\\)", out)))
  expect_true(any(grepl("Edges: 90", out)))
})

test_that("netcensus wraps a dynamic network and finds its window", {
  obs <- netcensus(nw_dyn)
  expect_true(obs$dynamic)
  expect_equal(obs$window, c(0, 30))
  expect_equal(obs$summary$edges.ever, network.edgecount(nw_dyn))
  ea <- obs$summary$edges.active
  expect_equal(names(ea), as.character(0:30))
  expect_equal(unname(ea["10"]), nrow(active_at(nw_dyn, 10)))
  out <- capture.output(print(obs))
  expect_true(any(grepl("dynamic \\(networkDynamic\\)", out)))
  expect_true(any(grepl("Observation window: 0 to 30", out)))

  # without net.obs.period, the window is the range of finite spell times
  nd <- nw_dyn
  network::delete.network.attribute(nd, "net.obs.period")
  obs2 <- netcensus(nd)
  sp <- networkDynamic::get.edge.activity(nd, as.spellList = TRUE)
  expect_equal(obs2$window[1], min(sp$onset[is.finite(sp$onset)]))
  expect_equal(obs2$window[2], max(sp$terminus[is.finite(sp$terminus)]))

  # an explicit window overrides
  obs3 <- netcensus(nw_dyn, window = c(5, 20))
  expect_equal(obs3$window, c(5, 20))
  expect_equal(names(obs3$summary$edges.active), as.character(5:20))
})

test_that("netcensus drops temporally extended vertex attributes", {
  nd <- nw_dyn
  nd <- networkDynamic::activate.vertex.attribute(nd, "status", "s",
                                                  onset = 0, terminus = Inf)
  expect_true("status.active" %in% list.vertex.attributes(nd))
  expect_message(obs <- netcensus(nd), "status.active")
  expect_false("status.active" %in% list.vertex.attributes(obs$newnetwork))
  expect_equal(obs$summary$attributes, "sex")
})

test_that("netcensus validates its inputs", {
  expect_error(netcensus(list()), "class `network`")
  expect_error(netcensus(network::network.initialize(10, directed = TRUE)),
               "undirected")
  expect_error(netcensus(network_initialize(10)), "no edges")
  expect_error(netcensus(nw_static, window = c(0, 10)), "networkDynamic")
  expect_error(netcensus(nw_dyn, window = c(10, 0)), "start < end")
  expect_error(netcensus(nw_dyn, window = 1), "start < end")
})

test_that("netdx refuses a netcensus object", {
  expect_error(netdx(netcensus(nw_static), nsims = 1, nsteps = 5),
               "does not apply")
})

test_that("make_stats_table accepts NULL targets", {
  stats <- list(matrix(c(1, 2, 3), ncol = 1, dimnames = list(NULL, "edges")))
  tab <- make_stats_table(stats, NULL)
  expect_true(is.na(tab$Target))
  expect_equal(tab[["Sim Mean"]], 2)
})


# Simulation: static census -----------------------------------------------------

test_that("netsim runs on a static census, tergmLite and networkDynamic", {
  obs <- netcensus(nw_static)
  param <- param.net(inf.prob = 0.3, act.rate = 1)
  for (tergmLite in c(TRUE, FALSE)) {
    control <- control.net(type = "SI", nsteps = 15, nsims = 1,
                           tergmLite = tergmLite,
                           resimulate.network = tergmLite, verbose = FALSE,
                           save.transmat = TRUE, save.run = TRUE)
    sim <- netsim(obs, param, init, control)
    expect_s3_class(sim, "netsim")
    expect_s3_class(sim$nwparam[[1]], "netcensus")
    test_net(sim)

    # the edge set is the observed one throughout
    ns <- get_nwstats(sim, network = 1)
    expect_equal(nrow(ns), 15)
    expect_true(all(ns$edges == 90))
    el <- if (tergmLite) sim$run[[1]]$el[[1]] else {
      active_at(get_network(sim), 15)
    }
    expect_equal(unclass(el)[, 1:2], unclass(as.edgelist(nw_static))[, 1:2],
                 ignore_attr = TRUE)

    # every transmission crossed an observed edge
    tm <- get_transmat(sim)
    obs_el <- as.edgelist(nw_static)
    keys <- paste(obs_el[, 1], obs_el[, 2])
    expect_true(all(paste(pmin(tm$sus, tm$inf), pmax(tm$sus, tm$inf)) %in% keys))

    # the nodal attribute came through
    expect_equal(sim$run[[1]]$attr$sex, rep(c("F", "M"), N / 2))

    # print and plot show a target of NA for the observed layer
    expect_output(print(sim), "netcensus")
    plot(sim)
    plot(sim, type = "formation")
  }
})


# Simulation: dynamic census ----------------------------------------------------

test_that("netsim reads the active edges of a dynamic census each step", {
  skip_on_cran()
  obs <- netcensus(nw_dyn)
  param <- param.net(inf.prob = 0.3, act.rate = 1)
  for (tergmLite in c(TRUE, FALSE)) {
    control <- control.net(type = "SI", nsteps = 25, nsims = 1,
                           tergmLite = tergmLite,
                           resimulate.network = tergmLite, verbose = FALSE,
                           save.transmat = TRUE, save.run = TRUE,
                           cumulative.edgelist = TRUE,
                           save.cumulative.edgelist = TRUE,
                           truncate.el.cuml = Inf)
    sim <- netsim(obs, param, init, control)
    test_net(sim)

    # the recorded edge count at each step is the census count at that time
    ns <- get_nwstats(sim, network = 1)
    census_counts <- vapply(1:25, function(t) nrow(active_at(nw_dyn, t)),
                            numeric(1))
    expect_equal(ns$edges, census_counts)

    # under tergmLite the final edgelist is the census at the last step
    if (tergmLite) {
      el <- sim$run[[1]]$el[[1]]
      expect_equal(unclass(el)[, 1:2], active_at(nw_dyn, 25),
                   ignore_attr = TRUE)
      expect_equal(attr(el, "n"), N)
    } else {
      # without tergmLite the observed spells are untouched
      expect_identical(networkDynamic::get.edge.activity(get_network(sim)),
                       networkDynamic::get.edge.activity(nw_dyn))
    }

    # every transmission crossed an edge active at that time
    tm <- get_transmat(sim)
    ok <- vapply(seq_len(nrow(tm)), function(i) {
      el <- active_at(nw_dyn, tm$at[i])
      any(el[, 1] == min(tm$sus[i], tm$inf[i]) &
            el[, 2] == max(tm$sus[i], tm$inf[i]))
    }, logical(1))
    expect_true(all(ok))

    # the cumulative edgelist follows the census: an edge that dissolved in
    # the census has a stop time
    cel <- sim$cumulative.edgelist[[1]]
    expect_true(any(!is.na(cel$stop)))
    expect_true(all(cel$start >= 0))
  }
})

test_that("a dynamic census works next to estimated and clique layers", {
  skip_on_cran()
  nw <- network_initialize(N)
  nw <- set_vertex_attribute(nw, "sex", rep(c("F", "M"), N / 2))
  nw <- set_vertex_attribute(nw, "hh_id", rep(1:20, each = 3))
  est <- netest(nw, ~edges, target.stats = 30,
                coef.diss = dissolution_coefs(~offset(edges), 10),
                verbose = FALSE)
  hh <- netclique(nw, group.attr = "hh_id")
  obs <- netcensus(nw_dyn)
  param <- param.net(inf.prob = multilayer(0.1, 0.3, 0.2), act.rate = 1)
  for (tergmLite in c(TRUE, FALSE)) {
    control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                           tergmLite = tergmLite, resimulate.network = TRUE,
                           verbose = FALSE, save.transmat = TRUE,
                           save.run = TRUE)
    sim <- netsim(list(est, obs, hh), param, init, control)
    test_net(sim)
    expect_s3_class(sim$nwparam[[2]], "netcensus")
    expect_s3_class(sim$nwparam[[3]], "netclique")
    # only the estimated layer has a coefficient
    expect_length(sim$nwparam[[1]]$coef.form, 1)
    expect_null(sim$nwparam[[2]]$coef.form)
    tm <- get_transmat(sim)
    expect_true(all(tm$network %in% 1:3))
    expect_true(all(tm$transProb[tm$network == 2] == 0.3))
    ns <- get_nwstats(sim, network = 2)
    expect_equal(ns$edges,
                 vapply(1:20, function(t) nrow(active_at(nw_dyn, t)),
                        numeric(1)))
  }
})

test_that("a census layer refuses vital dynamics and node changes", {
  obs <- netcensus(nw_dyn)
  param <- param.net(inf.prob = 0.3, a.rate = 0.01, ds.rate = 0.01,
                     di.rate = 0.01)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)
  expect_error(netsim(obs, param, init, control), "fixed node set")

  # a custom module that adds nodes is stopped in arrive_nodes
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  dat <- crosscheck.net(obs, param.net(inf.prob = 0.3), init, control)
  dat <- initialize.net(obs, param.net(inf.prob = 0.3), init, control, 1)
  expect_error(arrive_nodes(dat, 1), "fixed node set")
  expect_error(depart_nodes(dat, 1), "fixed node set")
  expect_s3_class(arrive_nodes(dat, 0), "netsim_dat")
})

test_that("duration tracking and runs past the window are caught", {
  obs <- netcensus(nw_dyn)
  param <- param.net(inf.prob = 0.3)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         tergmLite.track.duration = TRUE, verbose = FALSE)
  expect_error(netsim(obs, param, init, control), "track.duration")

  control <- control.net(type = "SI", nsteps = 40, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)
  expect_warning(sim <- netsim(obs, param, init, control),
                 "observation window")
  ns <- get_nwstats(sim, network = 1)
  # past the window the edge set is frozen
  expect_true(all(ns$edges[31:40] == ns$edges[31]))

  # a static census has no window and no warning
  control <- control.net(type = "SI", nsteps = 40, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)
  expect_silent(sim <- netsim(netcensus(nw_static), param, init, control))
})

test_that("a census layer survives a restart", {
  skip_on_cran()
  obs <- netcensus(nw_dyn)
  param <- param.net(inf.prob = 0.3)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE, save.run = TRUE)
  sim <- netsim(obs, param, init, control)
  control$start <- 11
  control$nsteps <- 20
  sim2 <- netsim(sim, param, init, control)
  expect_equal(nrow(sim2$epi$num), 20)
  ns <- get_nwstats(sim2, network = 1)
  expect_equal(ns$edges[11:20],
               vapply(11:20, function(t) nrow(active_at(nw_dyn, t)),
                      numeric(1)))
  test_net(sim2)
})
