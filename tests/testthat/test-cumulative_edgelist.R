context("Cumulative Edgelist (All SOC)")

nw <- network_initialize(n = 100)
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
  verbose = FALSE
)

param <- param.net(
  inf.prob = 0.3,
  act.rate = 0.1
)

init <- init.net(i.num = 10)

test_that("netsim, SI, Cumulative Edgelist", {
  skip_on_cran()
  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = NULL,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])

  expect_is(d, "data.frame")
  expect_equal(colnames(d), c("head", "tail", "start", "stop", "network"))

  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = 40,
    cumulative.edgelist = TRUE,
    save.cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])
  expect_gte(min(d$stop, na.rm = TRUE), 40)

  d <- get_partners(dat = mod[[1]], index_posit_ids = 1:10, networks = 1)
  expect_equal(colnames(d), c("index", "partner", "start", "stop", "network"))

  # expect errors when trying to access networks that don't exist
  expect_error(
    get_partners(dat = mod[[1]], index_posit_ids = 1:10, networks = 2)
  )

  expect_error(
    update_cumulative_edgelist(dat = mod[[1]], network = 2)
  )

  # Cumulative edgelist extraction and reachable set
  sim <- process_out.net(mod)
  el_cuml <- sim$cumulative.edgelist[[1]]
  el_cuml <- dedup_cumulative_edgelist(el_cuml)

  nnodes <- max(el_cuml$head, el_cuml$tail)

  from_step <- min(el_cuml$start)
  to_step <- 100
  nodes <- sample(nnodes, 10)
  # Only test that the functions run without error. See the functions examples
  # to check correctness
  el_tp <- get_forward_reachable(el_cuml, from_step, to_step, nodes)
  el_tp <- get_backward_reachable(el_cuml, from_step, to_step, nodes, "yes")
})

test_that("netsim, SI, Cumulative Edgelist - missing args", {
  skip_on_cran()
  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 10,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)

  expect_error(get_cumulative_edgelists_df(mod[[1]]))

  expect_equal(mod[[1]]$control$cumulative.edgelist, FALSE)
  expect_equal(mod[[1]]$control$truncate.el.cuml, 0)

})

test_that("netsim non-tergmLite cumulative edgelist matches networkDynamic (#1016)", {
  skip_on_cran()
  set.seed(1)
  nw <- network_initialize(n = 500)
  est <- netest(
    nw,
    formation = ~edges + concurrent,
    target.stats = c(0.75 * 500 / 2, 0.08 * 500),
    coef.diss = dissolution_coefs(~offset(edges), duration = 80),
    verbose = FALSE
  )

  control <- control.net(
    type = "SIS", nsims = 1, nsteps = 50,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = 99999,
    save.cumulative.edgelist = TRUE,
    save.network = TRUE,
    verbose = FALSE
  )

  set.seed(1)
  mod <- netsim(est, param.net(inf.prob = 0.5, act.rate = 2, rec.rate = 0.05),
                init.net(i.num = 25), control)

  cel_nd <- as.data.frame(get_network(mod, sim = 1))
  cel_epi <- mod$cumulative.edgelist$sim1

  # Row-count parity (the headline #1016 regression check).
  expect_equal(nrow(cel_epi), nrow(cel_nd))

  # Active-set parity at t=0 and t=1. networkDynamic spells are half-open
  # [onset, terminus); the cumulative edgelist is closed [start, stop] with
  # NA stop meaning ongoing. So an edge active at time t in nD has
  # onset <= t < terminus, and in the cumulative edgelist has
  # start <= t and (is.na(stop) | stop >= t).
  n_active_t0_nd  <- sum(cel_nd$onset <= 0 & cel_nd$terminus > 0)
  n_active_t0_epi <- sum(cel_epi$start <= 0 &
                           (is.na(cel_epi$stop) | cel_epi$stop >= 0))
  expect_equal(n_active_t0_epi, n_active_t0_nd)

  n_active_t1_nd  <- sum(cel_nd$onset <= 1 & cel_nd$terminus > 1)
  n_active_t1_epi <- sum(cel_epi$start <= 1 &
                           (is.na(cel_epi$stop) | cel_epi$stop >= 1))
  expect_equal(n_active_t1_epi, n_active_t1_nd)

  # Initial duration-1 edges (onset=0, terminus=1 in nD) -> start=0, stop=0.
  n_dur1_initial_nd  <- sum(cel_nd$onset == 0 & cel_nd$terminus == 1)
  n_dur1_initial_epi <- sum(cel_epi$start == 0 & cel_epi$stop == 0,
                            na.rm = TRUE)
  expect_equal(n_dur1_initial_epi, n_dur1_initial_nd)

  # Cross-section persistent edges live at start=0; nothing should be
  # earlier than that.
  expect_equal(min(cel_epi$start), 0)
})

test_that("truncate.el.cuml = 0 skips hist seed (#1016)", {
  skip_on_cran()
  set.seed(2)
  nw <- network_initialize(n = 100)
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 30,
    coef.diss = dissolution_coefs(~offset(edges), duration = 2),
    verbose = FALSE
  )

  control <- control.net(
    type = "SI", nsims = 1, nsteps = 3,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = 0,
    save.cumulative.edgelist = TRUE,
    verbose = FALSE
  )

  mod <- netsim(est, param.net(inf.prob = 0.1),
                init.net(i.num = 5), control)
  el_cuml <- mod$cumulative.edgelist$sim1

  # truncate=0 keeps only active edges: all rows must have stop==NA, and
  # in particular the start=0/stop=0 hist seed must not be applied.
  expect_true(all(is.na(el_cuml$stop)))
})

test_that("tergmLite cumulative edgelist captures initial-step edges (#1016)", {
  skip_on_cran()
  set.seed(3)
  # Short duration so a substantial fraction of cross-section edges dissolves
  # during the initial TERGM step. Without the fix in #1017 those rows would
  # be missing entirely; with option (b) they should be present as
  # start=0/stop=0 rows.
  nw <- network_initialize(n = 200)
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 80,
    coef.diss = dissolution_coefs(~offset(edges), duration = 2),
    verbose = FALSE
  )

  control <- control.net(
    type = "SI", nsims = 1, nsteps = 5,
    resimulate.network = TRUE, tergmLite = TRUE,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = Inf,
    save.cumulative.edgelist = TRUE,
    verbose = FALSE
  )

  set.seed(3)
  mod <- netsim(est, param.net(inf.prob = 0.1),
                init.net(i.num = 5), control)
  el_cuml <- mod$cumulative.edgelist$sim1

  # Cross-section persistent edges live at start=0 in tergmLite too, and
  # tergmLite is no longer limited on the initial-step dissolutions because
  # sim_nets_t1 stashes the t=0 edgelist before networkLite discards it.
  expect_equal(min(el_cuml$start), 0)
  expect_gt(sum(el_cuml$start == 0 & !is.na(el_cuml$stop) &
                  el_cuml$stop == 0), 0)
})

test_that("cumulative edgelist start values align with networkDynamic onsets (#1016)", {
  # Stronger pointwise check: every row's (start, stop) maps to a unique nD
  # spell (onset = start, terminus = stop + 1 with NA stop -> Inf).
  skip_on_cran()
  set.seed(4)
  nw <- network_initialize(n = 200)
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 60,
    coef.diss = dissolution_coefs(~offset(edges), duration = 5),
    verbose = FALSE
  )

  control <- control.net(
    type = "SI", nsims = 1, nsteps = 10,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = Inf,
    save.cumulative.edgelist = TRUE,
    save.network = TRUE,
    verbose = FALSE
  )

  set.seed(4)
  mod <- netsim(est, param.net(inf.prob = 0.1),
                init.net(i.num = 5), control)

  cel_nd  <- as.data.frame(get_network(mod, sim = 1))
  cel_epi <- mod$cumulative.edgelist$sim1

  # Build comparable spell signatures. nD's terminus is exclusive; ours is
  # inclusive, so terminus = stop + 1 with NA stop -> Inf. Right-censored
  # terminus values (edges ongoing at end of observation) likewise map to
  # Inf so they compare equal to NA stops in the cumulative edgelist.
  nd_term  <- ifelse(cel_nd$terminus.censored, Inf, cel_nd$terminus)
  epi_term <- ifelse(is.na(cel_epi$stop), Inf, cel_epi$stop + 1)
  nd_key <- sort(paste(pmin(cel_nd$tail, cel_nd$head),
                       pmax(cel_nd$tail, cel_nd$head),
                       cel_nd$onset, nd_term, sep = "_"))
  epi_key <- sort(paste(pmin(cel_epi$head, cel_epi$tail),
                        pmax(cel_epi$head, cel_epi$tail),
                        cel_epi$start, epi_term, sep = "_"))
  expect_identical(epi_key, nd_key)
})

test_that("netsim, SI, Cumulative Edgelist with arrivals and departures", {
  skip_on_cran()
  nw <- network_initialize(n = 100)

  formation <- ~edges
  target.stats <- 70
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38, d.rate = 0.01)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, act.rate = 1,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)

  init <- init.net(i.num = 10)
  control <- control.net(
    type = "SI", nsteps = 100, nsims = 1,
    tergmLite = TRUE, # required as we want the removal of inactive nodes
    resimulate.network = TRUE, verbose = FALSE,
    truncate.el.cuml = Inf,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  sim1 <- netsim(est1, param, init, control)

  dat <- sim1[[1]]
  uid <- get_unique_ids(dat)
  pid <- get_posit_ids(dat)
  is_active_posit_ids(dat, pid)

  el_cuml <- get_cumulative_edgelists_df(dat)
  spids <- sample(pid, 20)
  partners <- get_partners(dat, spids)
  degree_cuml <- get_cumulative_degree(dat, spids)

  # checks the unicity of UIDs
  expect_true(all(table(uid) == 1))
  # checks that all attributed UID are not in the currently attributed UIDs
  expect_false(setequal(seq_len(dat$run$last_unique_id), get_unique_ids(dat)))
  # checks that all attributed UID are not in the currently attributed PIDs
  expect_false(setequal(seq_len(dat$run$last_unique_id), get_posit_ids(dat)))

  # checks the translation between UID and PID in get_partners
  expect_true(all(get_posit_ids(dat, unique(partners$index)) %in% spids))

})

test_that("dedup_cumulative_edgelist merges overlapping spells", {
  spells <- function(start, stop, head = 1, tail = 2) {
    data.frame(head = head, tail = tail, start = start, stop = stop)
  }

  ## overlapping spells of one edge collapse into a single row
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 2), c(3, 4))),
               spells(1, 4))
  ## a chain of overlaps collapses transitively
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 2, 4), c(3, 5, 8))),
               spells(1, 8))
  ## a spell fully contained in another is absorbed
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 2), c(9, 3))),
               spells(1, 9))

  ## disjoint spells are kept, including ones that merely touch
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 5), c(2, 6))),
               spells(c(1, 5), c(2, 6)))
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 3), c(2, 5))),
               spells(c(1, 3), c(2, 5)))

  ## the input need not be sorted
  expect_equal(dedup_cumulative_edgelist(spells(c(7, 1), c(9, 8))),
               spells(1, 9))

  ## edges are merged separately, and the output is sorted by edge then start
  el <- data.frame(head = c(1, 1, 1, 2), tail = c(2, 2, 2, 3),
                   start = c(2, 8, 1, 4), stop = c(5, 9, 3, 6))
  expect_equal(dedup_cumulative_edgelist(el),
               data.frame(head = c(1, 1, 2), tail = c(2, 2, 3),
                          start = c(1, 8, 4), stop = c(5, 9, 6)))
})

test_that("dedup_cumulative_edgelist keeps open-ended spells", {
  spells <- function(start, stop, head = 1, tail = 2) {
    data.frame(head = head, tail = tail, start = start, stop = stop)
  }

  ## an NA stop means the edge never dissolved, so it is never dropped
  expect_equal(dedup_cumulative_edgelist(spells(1, NA_real_)), spells(1, NA_real_))

  ## an open-ended spell overlapping a closed one absorbs it and stays open
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 2), c(3, NA))),
               spells(1, NA_real_))
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 3), c(NA, 5))),
               spells(1, NA_real_))

  ## it also absorbs every later spell of that edge
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 3, 20), c(NA, 5, 25))),
               spells(1, NA_real_))
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 4), c(NA, NA))),
               spells(1, NA_real_))

  ## an open-ended spell starting after a closed one ends stays separate
  expect_equal(dedup_cumulative_edgelist(spells(c(1, 5), c(2, NA))),
               spells(c(1, 5), c(2, NA)))
})

test_that("dedup_cumulative_edgelist handles other columns and edge cases", {
  ## spells are kept apart by any other column, such as `network`
  el <- data.frame(head = c(1, 1), tail = c(2, 2), start = c(1, 2),
                   stop = c(3, 4), network = c(1, 2))
  expect_equal(dedup_cumulative_edgelist(el), el)
  el$network <- c(1, 1)
  expect_equal(dedup_cumulative_edgelist(el),
               data.frame(head = 1, tail = 2, start = 1, stop = 4,
                          network = 1))

  ## empty input, storage modes, and the object class are preserved
  empty <- data.frame(head = numeric(0), tail = numeric(0),
                      start = numeric(0), stop = numeric(0))
  expect_equal(dedup_cumulative_edgelist(empty), empty)

  int.el <- data.frame(head = c(1L, 1L), tail = c(2L, 2L),
                       start = c(1L, 2L), stop = c(3L, 4L))
  expect_type(dedup_cumulative_edgelist(int.el)$stop, "integer")

  class(int.el) <- c("cumulative_edgelist", class(int.el))
  expect_s3_class(dedup_cumulative_edgelist(int.el), "cumulative_edgelist")

  expect_error(
    dedup_cumulative_edgelist(data.frame(head = 1, tail = 2, start = NA,
                                         stop = 3)),
    "cannot contain NA"
  )
})
