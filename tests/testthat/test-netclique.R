context("Clique Network Layers")

# Shared fixtures ------------------------------------------------------------

# A small population built from household types, with an age group per person
# and a household id per household. Every child lives with at least one adult.
hh_types <- c("adult"                     = 0.15,
              "adult adult"               = 0.20,
              "elderly"                   = 0.12,
              "elderly elderly"           = 0.10,
              "adult elderly"             = 0.03,
              "adult child"               = 0.05,
              "adult adult child"         = 0.15,
              "adult adult child child"   = 0.15,
              "adult adult child elderly" = 0.05)

set.seed(11)
N <- 200
pop <- sample_groups(N, hh_types, attr.name = "age")
nw <- network_initialize(n = N)
nw <- set_vertex_attribute(nw, "age", pop$age)
nw <- set_vertex_attribute(nw, "hh_id", pop$group)

est_hh <- netclique(nw, group.attr = "hh_id", arrivals = "join")
est_com <- netest(nw, formation = ~edges + nodematch("age"),
                  target.stats = c(100, 70),
                  coef.diss = dissolution_coefs(~offset(edges), 20),
                  verbose = FALSE)
init <- init.net(i.num = 10)

# The clique invariant: every clique edge joins two members of the same
# group, and every pair of active co-members is an edge.
expect_clique_layer <- function(el, group, active = rep(1, length(group))) {
  in_group <- which(active == 1 & !is.na(group))
  sizes <- table(group[in_group])
  expect_equal(nrow(el), sum(choose(sizes, 2)))
  expect_true(all(group[el[, 1]] == group[el[, 2]]))
  expect_true(all(active[el[, 1]] == 1 & active[el[, 2]] == 1))
}


# Constructor ------------------------------------------------------------------

test_that("netclique builds cliques from a grouping attribute", {
  expect_s3_class(est_hh, "netclique")
  sizes <- tabulate(pop$group)
  expect_equal(est_hh$summary$edges, sum(choose(sizes, 2)))
  expect_equal(network.edgecount(est_hh$newnetwork), sum(choose(sizes, 2)))
  expect_equal(est_hh$summary$n, N)
  expect_equal(est_hh$summary$isolates, sum(sizes == 1))
  expect_equal(est_hh$target.stats, est_hh$summary$edges)
  expect_equal(est_hh$target.stats.names, "edges")
  expect_false(est_hh$edapprox)
  expect_null(est_hh$formation)
  expect_null(est_hh$coef.diss)
  expect_equal(est_hh$group.attr, "hh_id")
  expect_equal(est_hh$arrivals, "join")

  el <- as.edgelist(est_hh$newnetwork)
  expect_clique_layer(el, pop$group)

  # the layer network carries the nodal attributes of nw
  expect_equal(get_vertex_attribute(est_hh$newnetwork, "age"), pop$age)
  expect_equal(get_vertex_attribute(est_hh$newnetwork, "hh_id"), pop$group)
})

test_that("netclique treats NA group values as isolates", {
  g <- pop$group
  g[1:10] <- NA
  nw_na <- set_vertex_attribute(nw, "hh_id", g)
  es <- netclique(nw_na, group.attr = "hh_id")
  el <- as.edgelist(es$newnetwork)
  expect_false(any(el %in% 1:10))
  expect_clique_layer(el, g)
})

test_that("netclique accepts character group ids", {
  nw_chr <- set_vertex_attribute(nw, "hh_id", paste0("hh", pop$group))
  es <- netclique(nw_chr, group.attr = "hh_id")
  expect_equal(es$summary$edges, est_hh$summary$edges)
})

test_that("netclique validates its inputs", {
  expect_error(netclique(list(), group.attr = "hh_id"), "class `network`")
  expect_error(netclique(network::network.initialize(10, directed = TRUE),
                         group.attr = "hh_id"), "undirected")
  expect_error(netclique(nw), "single attribute")
  expect_error(netclique(nw, group.attr = c("a", "b")), "single attribute")
  expect_error(netclique(nw, group.attr = "nope"), "no vertex attribute")
  nw_lgl <- set_vertex_attribute(nw, "hh_id", rep(TRUE, N))
  expect_error(netclique(nw_lgl, group.attr = "hh_id"), "numeric or character")
  expect_error(netclique(nw, group.attr = "hh_id", arrivals.FUN = identity),
               "applies only with")
  expect_error(netclique(nw, group.attr = "hh_id", arrivals = "join",
                         arrivals.FUN = 1), "must be a function")
  # edges already on nw are ignored
  nw_e <- add.edges(nw, tail = c(1, 3), head = c(2, 4))
  expect_equal(netclique(nw_e, group.attr = "hh_id")$summary$edges,
               est_hh$summary$edges)
})

test_that("print.netclique reports the layer and mean degree by attribute", {
  out <- capture.output(print(est_hh))
  expect_true(any(grepl("Grouping attribute: hh_id", out)))
  expect_true(any(grepl("Group Size Distribution", out)))
  expect_true(any(grepl("Arrivals: join", out)))
  out <- capture.output(print(est_hh, by = "age"))
  expect_true(any(grepl("Mean Degree by `age`", out)))
  expect_error(print(est_hh, by = "nope"), "no vertex attribute")
})

test_that("netdx refuses a netclique object", {
  expect_error(netdx(est_hh, nsims = 1, nsteps = 5), "does not apply")
})


# Group construction helpers ----------------------------------------------------

test_that("sample_groups builds a population of exactly n from group types", {
  set.seed(2)
  p <- sample_groups(1000, hh_types, attr.name = "age")
  expect_equal(nrow(p), 1000)
  expect_equal(names(p), c("group", "type", "age"))
  # groups are consecutive and each is one drawn type
  expect_true(all(diff(p$group) %in% c(0, 1)))
  expect_equal(p$group[1], 1)
  members <- strsplit(names(hh_types), " ")
  by_group <- split(p, p$group)
  expect_true(all(vapply(by_group, function(g) {
    identical(g$age, members[[g$type[1]]])
  }, logical(1))))
  # every child shares a household with an adult, by construction of the types
  has_adult <- tapply(p$age == "adult", p$group, any)
  expect_true(all(has_adult[as.character(p$group[p$age == "child"])]))
  # the realized type frequencies follow the weights
  freq <- tabulate(p$type[!duplicated(p$group)], nbins = length(hh_types))
  expect_gt(cor(freq, hh_types), 0.9)
})

test_that("sample_groups accepts templates with prob, and data.frame types", {
  set.seed(3)
  p <- sample_groups(50, c("a", "a b"), prob = c(1, 1), attr.name = "role")
  expect_equal(nrow(p), 50)
  expect_true(all(p$role %in% c("a", "b")))

  types <- list(
    data.frame(age = "adult", sex = "F"),
    data.frame(age = c("adult", "adult"), sex = c("F", "M")),
    data.frame(age = c("adult", "adult", "child"), sex = c("F", "M", "F"))
  )
  p2 <- sample_groups(100, types, prob = c(0.3, 0.4, 0.3))
  expect_equal(nrow(p2), 100)
  expect_equal(names(p2), c("group", "type", "age", "sex"))
  expect_true(all(p2$sex[p2$age == "child"] == "F"))

  # the last group is truncated only when no type fits
  expect_message(p3 <- sample_groups(7, c("a a", "a a a a")), "truncated")
  expect_equal(nrow(p3), 7)
  expect_silent(p4 <- sample_groups(7, c("a", "a a a a")))
  expect_equal(nrow(p4), 7)
})

test_that("sample_groups validates its inputs", {
  expect_error(sample_groups(0, "a"), "positive integer")
  expect_error(sample_groups(10, c(a = "a")), "named character")
  expect_error(sample_groups(10, c(0.5, 0.5)), "must be named")
  expect_error(sample_groups(10, c(a = 0.5), prob = 1), "not both")
  expect_error(sample_groups(10, "a", prob = c(1, 1)), "one non-negative weight")
  expect_error(sample_groups(10, list(1)), "one `data.frame`")
  expect_error(sample_groups(10, list(data.frame(group = 1))), "cannot be named")
  expect_error(sample_groups(10, 1:2), "must be named")
})

test_that("assign_groups keeps every dependent with an anchor", {
  set.seed(4)
  n <- 2000
  age <- sample(c("child", "adult", "elderly"), n, replace = TRUE,
                prob = c(0.22, 0.60, 0.18))
  size.dist <- c("1" = 0.28, "2" = 0.34, "3" = 0.15, "4" = 0.14, "5" = 0.09)
  g <- assign_groups(size.dist, role = age,
                     anchor = c("adult", "elderly"), dependent = "child")
  expect_length(g, n)
  expect_false(anyNA(g))
  expect_equal(sort(unique(g)), seq_len(max(g)))
  has_anchor <- tapply(age %in% c("adult", "elderly"), g, any)
  expect_true(all(has_anchor[as.character(g[age == "child"])]))
  # the realized size distribution is close to the target
  realized <- prop.table(table(factor(tabulate(g), levels = 1:5)))
  expect_lt(max(abs(as.numeric(realized) - size.dist)), 0.05)

  # the same population, with only adults as anchors
  g2 <- assign_groups(size.dist, role = age, anchor = "adult",
                      dependent = "child")
  has_adult <- tapply(age == "adult", g2, any)
  expect_true(all(has_adult[as.character(g2[age == "child"])]))
})

test_that("assign_groups without roles assigns random groups of the drawn sizes", {
  set.seed(5)
  g <- assign_groups(c(0.5, 0.5), n = 500)
  expect_length(g, 500)
  expect_true(all(tabulate(g) %in% 1:2))
  g2 <- assign_groups(c("3" = 1), n = 30)
  expect_true(all(tabulate(g2) == 3))
})

test_that("assign_groups handles too many dependents and no anchors", {
  set.seed(6)
  role <- rep(c("child", "adult"), c(90, 10))
  expect_message(g <- assign_groups(c("2" = 1), role = role, anchor = "adult",
                                    dependent = "child"),
                 "exceed the open slots")
  has_adult <- tapply(role == "adult", g, any)
  expect_true(all(has_adult[as.character(g[role == "child"])]))

  expect_message(g2 <- assign_groups(c("2" = 1), role = rep("child", 20),
                                     anchor = "adult", dependent = "child"),
                 "No anchors")
  expect_length(g2, 20)
})

test_that("assign_groups validates its inputs", {
  expect_error(assign_groups(c(-1, 1), n = 10), "non-negative")
  expect_error(assign_groups(c(a = 1), n = 10), "positive integer group sizes")
  expect_error(assign_groups(c(1, 1)), "population size")
  expect_error(assign_groups(c(1, 1), role = 1:5, n = 4), "does not match")
  expect_error(assign_groups(c(1, 1), role = 1:5, dependent = 1), "requires")
  expect_error(assign_groups(c(1, 1), role = 1:5, anchor = 1, dependent = 1),
               "must not share")
})

test_that("the two helpers compose into a netclique layer", {
  set.seed(7)
  age <- sample(c("child", "adult", "elderly"), 300, replace = TRUE,
                prob = c(0.25, 0.55, 0.20))
  hh <- assign_groups(c("1" = 0.3, "2" = 0.35, "3" = 0.2, "4" = 0.15),
                      role = age, anchor = c("adult", "elderly"),
                      dependent = "child")
  nw2 <- network_initialize(300)
  nw2 <- set_vertex_attribute(nw2, "age", age)
  nw2 <- set_vertex_attribute(nw2, "hh_id", hh)
  es <- netclique(nw2, group.attr = "hh_id")
  expect_clique_layer(as.edgelist(es$newnetwork), hh)
  # a child's household degree is at least one (the anchor)
  deg <- get_degree(as.edgelist(es$newnetwork))
  expect_true(all(deg[age == "child"] >= 1))
})


# Simulation: closed population -------------------------------------------------

test_that("netsim runs a clique plus TERGM model, tergmLite, with layer params", {
  param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE, save.run = TRUE)
  sim <- netsim(list(est_hh, est_com), param, init, control)
  expect_s3_class(sim, "netsim")
  expect_equal(sim$num.nw, 2)
  test_net(sim)

  # the clique layer's edgelist is unchanged after the run
  expect_equal(unclass(sim$run[[1]]$el[[1]])[, 1:2],
               unclass(as.edgelist(est_hh$newnetwork))[, 1:2],
               ignore_attr = TRUE)
  # the TERGM layer's coefficient is untouched by the edges correction in a
  # closed population, and the clique layer has none
  expect_equal(sim$nwparam[[2]]$coef.form, est_com$coef.form)
  expect_null(sim$nwparam[[1]]$coef.form)
  expect_s3_class(sim$nwparam[[1]], "netclique")

  # network statistics: edges of the static layer are constant at its count
  ns <- get_nwstats(sim, network = 1)
  expect_true(all(ns$edges == est_hh$summary$edges))
  ns2 <- get_nwstats(sim, network = 2)
  expect_equal(names(ns2), c("time", "sim", "edges", "nodematch.age"))

  # transmissions record their layer, and both layers transmit
  tm <- get_transmat(sim)
  expect_true("network" %in% names(tm))
  expect_true(all(tm$network %in% 1:2))
  expect_true(all(tm$transProb[tm$network == 1] == 0.3))
  expect_true(all(tm$transProb[tm$network == 2] == 0.05))

  # print and plot work for both layers
  expect_output(print(sim, network = 1), "Target")
  expect_output(print(sim, network = 2), "nodematch.age")
  plot(sim)
  plot(sim, type = "formation", network = 1)
  plot(sim, type = "formation", network = 2)
})

test_that("a clique layer may come after the TERGM layer, networkDynamic, SIS", {
  skip_on_cran()
  param <- param.net(inf.prob = multilayer(0.05, 0.3),
                     act.rate = multilayer(1, 2), rec.rate = 0.05)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1,
                         tergmLite = FALSE, resimulate.network = TRUE,
                         verbose = FALSE, save.diss.stats = TRUE,
                         save.transmat = TRUE)
  sim <- netsim(list(est_com, est_hh), param, init, control)
  expect_s3_class(sim, "netsim")
  test_net(sim)

  nd <- get_network(sim, network = 2)
  expect_s3_class(nd, "networkDynamic")
  el <- networkDynamic::get.dyads.active(nd, at = 10)
  expect_equal(nrow(el), est_hh$summary$edges)

  tm <- get_transmat(sim)
  expect_true(all(tm$actRate[tm$network == 2] == 2))
  expect_true(all(tm$actRate[tm$network == 1] == 1))

  # duration and dissolution diagnostics exist for the TERGM layer only
  expect_output(print(sim, network = 1), "Duration Statistics")
  expect_output(print(sim, network = 2), "netclique")
  plot(sim, type = "duration", network = 1)
  expect_error(plot(sim, type = "duration", network = 2), "netclique")
})

test_that("a clique layer works without network resimulation", {
  skip_on_cran()
  param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
  for (tergmLite in c(FALSE, TRUE)) {
    control <- suppressWarnings(
      control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = tergmLite,
                  resimulate.network = FALSE, verbose = FALSE)
    )
    sim <- netsim(list(est_hh, est_com), param, init, control)
    expect_s3_class(sim, "netsim")
    ns <- get_nwstats(sim, network = 1)
    expect_equal(nrow(ns), 10)
    expect_true(all(ns$edges == est_hh$summary$edges))
  }
})

test_that("a model may consist of clique layers only", {
  skip_on_cran()
  param <- param.net(inf.prob = 0.2, act.rate = 1)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  sim <- netsim(est_hh, param, init, control)
  expect_s3_class(sim, "netsim")
  expect_equal(sim$num.nw, 1)
  test_net(sim)

  # infections stay within households: every infected node either was a seed
  # or shares a household with an earlier infected node
  control <- control.net(type = "SI", nsteps = 30, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE)
  sim <- netsim(est_hh, param.net(inf.prob = 0.5), init, control)
  tm <- get_transmat(sim)
  expect_true(all(pop$group[tm$sus] == pop$group[tm$inf]))
})

test_that("a two-group model accepts multilayer inf.prob and inf.prob.g2", {
  skip_on_cran()
  nw2 <- set_vertex_attribute(nw, "group", rep(1:2, length.out = N))
  es <- netclique(nw2, group.attr = "hh_id")
  ec <- netest(nw2, formation = ~edges, target.stats = 100,
               coef.diss = dissolution_coefs(~offset(edges), 20),
               verbose = FALSE)
  param <- param.net(inf.prob = multilayer(0.3, 0.05),
                     inf.prob.g2 = multilayer(0.2, 0.02), act.rate = 1)
  init2 <- init.net(i.num = 5, i.num.g2 = 5)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE)
  sim <- netsim(list(es, ec), param, init2, control)
  expect_s3_class(sim, "netsim")
  test_net(sim)
  tm <- get_transmat(sim)
  grp <- rep(1:2, length.out = N)
  expect_true(all(tm$transProb[tm$network == 1 & grp[tm$sus] == 1] == 0.3))
  expect_true(all(tm$transProb[tm$network == 1 & grp[tm$sus] == 2] == 0.2))
  expect_true(all(tm$transProb[tm$network == 2 & grp[tm$sus] == 1] == 0.05))
  expect_true(all(tm$transProb[tm$network == 2 & grp[tm$sus] == 2] == 0.02))
})

test_that("multilayer parameters may also be time-varying", {
  skip_on_cran()
  param <- param.net(inf.prob = multilayer(c(0.5, 0.1), 0.05), act.rate = 1)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE)
  sim <- netsim(list(est_hh, est_com), param, init, control)
  tm <- get_transmat(sim)
  hh_tm <- tm[tm$network == 1, ]
  expect_true(all(hh_tm$transProb[hh_tm$infDur <= 1] == 0.5))
  expect_true(all(hh_tm$transProb[hh_tm$infDur > 1] == 0.1))
})

test_that("crosscheck.net rejects mismatched layers and parameters", {
  param <- param.net(inf.prob = multilayer(0.3, 0.05, 0.1), act.rate = 1)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  expect_error(netsim(list(est_hh, est_com), param, init, control),
               "multilayer parameter `inf.prob` has length 3")

  param <- param.net(inf.prob = 0.3, act.rate = 1)
  nw_small <- set_vertex_attribute(network_initialize(10), "hh_id", rep(1:5, 2))
  es_small <- netclique(nw_small, group.attr = "hh_id")
  expect_error(netsim(list(es_small, est_com), param, init, control),
               "same number of nodes")

  expect_error(netsim(list(est_hh, "not a layer"), param, init, control),
               "netest or netclique")
})

test_that("print.param.net shows multilayer parameters", {
  param <- param.net(inf.prob = multilayer(0.3, 0.05),
                     act.rate = multilayer(c(1, 2), 1))
  out <- capture.output(print(param))
  expect_true(any(grepl("inf.prob = multilayer(0.3, 0.05)", out, fixed = TRUE)))
  expect_true(any(grepl("act.rate = multilayer(c(1, 2), 1)", out, fixed = TRUE)))
})


# Simulation: open population and arrival rules ----------------------------------

param_open <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1,
                        a.rate = 0.05, ds.rate = 0.02, di.rate = 0.02)

test_that("join: arrivals join a group in proportion to its size, tergmLite", {
  skip_on_cran()
  set.seed(21)
  control <- control.net(type = "SI", nsteps = 25, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE, cumulative.edgelist = TRUE,
                         save.cumulative.edgelist = TRUE,
                         truncate.el.cuml = Inf,
                         tergmLite.track.duration = TRUE)
  sim <- netsim(list(est_hh, est_com), param_open, init, control)
  test_net(sim)
  expect_gt(sum(sim$epi$a.flow[, 1], na.rm = TRUE), 0)
  expect_gt(sum(sim$epi$ds.flow[, 1] + sim$epi$di.flow[, 1], na.rm = TRUE), 0)

  run <- sim$run[[1]]
  hh <- run$attr$hh_id
  el <- run$el[[1]]
  expect_equal(attr(el, "n"), length(hh))
  expect_clique_layer(el, hh)

  # every arrival was placed in a group that existed before it
  arrived <- which(run$attr$entrTime > 1)
  expect_gt(length(arrived), 0)
  expect_false(anyNA(hh[arrived]))
  expect_true(all(hh[arrived] %in% pop$group))

  # the TERGM layer received the edges correction, the clique layer did not
  expect_equal(sim$nwparam[[2]]$coef.form[1],
               est_com$coef.form[1] + log(N) - log(run$num),
               tolerance = 1e-6)
  expect_null(sim$nwparam[[1]]$coef.form)

  # duration tracking: original edges toggled at 0, arrival edges later
  lt <- run$net_attr[[1]]$lasttoggle
  expect_equal(nrow(lt), nrow(el))
  expect_true(all(lt[, 3] >= 0 & lt[, 3] <= 25))
  expect_equal(run$net_attr[[1]]$time, 25)

  # cumulative edgelist: the original clique edges start at 0
  cel <- sim$cumulative.edgelist[[1]]
  cel1 <- cel[cel$network == 1, ]
  expect_true(all(cel1$start >= 0))
  expect_equal(sum(cel1$start == 0), est_hh$summary$edges)
})

test_that("join: arrivals join a group, networkDynamic", {
  skip_on_cran()
  set.seed(22)
  control <- control.net(type = "SI", nsteps = 15, nsims = 1, tergmLite = FALSE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  sim <- netsim(list(est_hh, est_com), param_open, init, control)
  test_net(sim)
  run <- sim$run[[1]]
  hh <- run$attr$hh_id
  active <- run$attr$active
  arrived <- which(run$attr$entrTime > 1)
  expect_gt(length(arrived), 0)
  expect_false(anyNA(hh[arrived]))

  nd <- get_network(sim, network = 1)
  el <- networkDynamic::get.dyads.active(nd, at = 15)
  expect_clique_layer(el, hh, active)
})

test_that("new: each arrival starts a group of its own", {
  skip_on_cran()
  set.seed(23)
  es <- netclique(nw, group.attr = "hh_id", arrivals = "new")
  for (tergmLite in c(TRUE, FALSE)) {
    control <- control.net(type = "SI", nsteps = 15, nsims = 1,
                           tergmLite = tergmLite, resimulate.network = TRUE,
                           verbose = FALSE, save.run = TRUE)
    sim <- netsim(list(es, est_com), param_open, init, control)
    run <- sim$run[[1]]
    hh <- run$attr$hh_id
    arrived <- which(run$attr$entrTime > 1)
    expect_gt(length(arrived), 0)
    expect_false(anyNA(hh[arrived]))
    expect_true(all(hh[arrived] > max(pop$group)))
    expect_false(any(duplicated(hh[arrived])))
    el <- if (tergmLite) run$el[[1]] else {
      networkDynamic::get.dyads.active(get_network(sim, network = 1), at = 15)
    }
    expect_false(any(arrived %in% el))
    expect_clique_layer(el, hh, run$attr$active)
  }
})

test_that("new: character group ids get unique new values", {
  skip_on_cran()
  set.seed(24)
  nw_chr <- set_vertex_attribute(nw, "hh_id", paste0("hh", pop$group))
  es <- netclique(nw_chr, group.attr = "hh_id", arrivals = "new")
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  sim <- netsim(list(es, est_com), param_open, init, control)
  hh <- sim$run[[1]]$attr$hh_id
  arrived <- which(sim$run[[1]]$attr$entrTime > 1)
  expect_gt(length(arrived), 0)
  expect_true(is.character(hh))
  expect_false(any(hh[arrived] %in% paste0("hh", pop$group)))
  expect_false(any(duplicated(hh[arrived])))
})

test_that("isolate: arrivals get no edges and an NA group", {
  skip_on_cran()
  set.seed(25)
  es <- netclique(nw, group.attr = "hh_id", arrivals = "isolate")
  control <- control.net(type = "SI", nsteps = 15, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  sim <- netsim(list(es, est_com), param_open, init, control)
  run <- sim$run[[1]]
  arrived <- which(run$attr$entrTime > 1)
  expect_gt(length(arrived), 0)
  expect_true(all(is.na(run$attr$hh_id[arrived])))
  expect_false(any(arrived %in% run$el[[1]]))
  expect_clique_layer(run$el[[1]], run$attr$hh_id)

  # a user attr.rules entry for the group attribute is respected
  control$attr.rules <- list(hh_id = 1)
  sim <- netsim(list(es, est_com), param_open, init, control)
  run <- sim$run[[1]]
  arrived <- which(run$attr$entrTime > 1)
  expect_true(all(run$attr$hh_id[arrived] == 1))
})

test_that("join with arrivals.FUN: the user function picks the group", {
  skip_on_cran()
  set.seed(26)
  # newborns go to a household that already has a child
  es <- netclique(nw, group.attr = "hh_id", arrivals = "join",
    arrivals.FUN = function(dat, at, new_ids, network) {
      hh <- get_attr(dat, "hh_id")
      age <- get_attr(dat, "age")
      entr <- get_attr(dat, "entrTime")
      pool <- hh[which(age == "child" & !is.na(hh) & entr < at)]
      pool[sample.int(length(pool), length(new_ids), replace = TRUE)]
    })
  control <- control.net(type = "SI", nsteps = 15, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE, attr.rules = list(age = "child"))
  sim <- netsim(list(es, est_com), param_open, init, control)
  run <- sim$run[[1]]
  hh <- run$attr$hh_id
  arrived <- which(run$attr$entrTime > 1)
  expect_gt(length(arrived), 0)
  expect_false(anyNA(hh[arrived]))
  expect_clique_layer(run$el[[1]], hh)
  # every household that received an arrival had a child in it originally
  hh_with_child <- unique(pop$group[pop$age == "child"])
  expect_true(all(hh[arrived] %in% hh_with_child))

  # an NA from the function leaves the node isolated
  es_na <- netclique(nw, group.attr = "hh_id", arrivals = "join",
    arrivals.FUN = function(dat, at, new_ids, network) {
      rep(NA, length(new_ids))
    })
  sim <- netsim(list(es_na, est_com), param_open, init, control)
  run <- sim$run[[1]]
  arrived <- which(run$attr$entrTime > 1)
  expect_true(all(is.na(run$attr$hh_id[arrived])))
  expect_false(any(arrived %in% run$el[[1]]))

  # a function returning the wrong length is an error
  es_bad <- netclique(nw, group.attr = "hh_id", arrivals = "join",
    arrivals.FUN = function(dat, at, new_ids, network) 1)
  expect_error(netsim(list(es_bad, est_com), param_open, init, control),
               "returned 1 group ids")
})

test_that("two arrivals joining the same group in one step are connected", {
  skip_on_cran()
  # force many arrivals into a single group every step
  es <- netclique(nw, group.attr = "hh_id", arrivals = "join",
    arrivals.FUN = function(dat, at, new_ids, network) {
      rep(1, length(new_ids))
    })
  param <- param.net(inf.prob = 0.1, act.rate = 1, a.rate = 0.1,
                     ds.rate = 0.001, di.rate = 0.001)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  set.seed(27)
  sim <- netsim(list(es, est_com), param, init, control)
  run <- sim$run[[1]]
  hh <- run$attr$hh_id
  arrived <- which(run$attr$entrTime > 1)
  expect_gt(length(arrived), 1)
  expect_true(all(hh[arrived] == 1))
  members <- which(hh == 1)
  el <- run$el[[1]]
  in_g1 <- el[hh[el[, 1]] == 1, , drop = FALSE]
  expect_equal(nrow(in_g1), choose(length(members), 2))
  expect_clique_layer(el, hh)
})

test_that("the group attribute is copied from the clique layer when absent", {
  skip_on_cran()
  # the TERGM layer's network lacks hh_id; netsim takes it from the clique layer
  nw_no_hh <- network_initialize(N)
  nw_no_hh <- set_vertex_attribute(nw_no_hh, "age", pop$age)
  ec <- netest(nw_no_hh, formation = ~edges, target.stats = 100,
               coef.diss = dissolution_coefs(~offset(edges), 20),
               verbose = FALSE)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  set.seed(28)
  sim <- netsim(list(ec, est_hh), param_open, init, control)
  run <- sim$run[[1]]
  expect_true("hh_id" %in% names(run$attr))
  expect_clique_layer(run$el[[2]], run$attr$hh_id)
})

test_that("a clique layer survives a restart", {
  skip_on_cran()
  set.seed(29)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE)
  sim <- netsim(list(est_hh, est_com), param_open, init, control)
  control$start <- 11
  control$nsteps <- 20
  sim2 <- netsim(sim, param_open, init, control)
  expect_equal(nrow(sim2$epi$num), 20)
  expect_s3_class(sim2$nwparam[[1]], "netclique")
  run <- sim2$run[[1]]
  expect_clique_layer(run$el[[1]], run$attr$hh_id)
  test_net(sim2)
})

test_that("multi-layer transmissions are not attributed to the first layer only", {
  skip_on_cran()
  # two identical clique layers with equal transmission probability: a node
  # exposed on both is credited to either with equal chance, so the recorded
  # shares are close to even rather than all on layer 1
  es1 <- netclique(nw, group.attr = "hh_id")
  param <- param.net(inf.prob = multilayer(0.4, 0.4), act.rate = 1)
  control <- control.net(type = "SI", nsteps = 30, nsims = 5, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE)
  set.seed(30)
  sim <- netsim(list(es1, es1), param, init, control)
  tm <- do.call(rbind, lapply(1:5, function(s) get_transmat(sim, sim = s)))
  share <- mean(tm$network == 1)
  expect_gt(share, 0.35)
  expect_lt(share, 0.65)
})

test_that("netsim objects with a clique layer merge and truncate", {
  skip_on_cran()
  param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.transmat = TRUE)
  set.seed(31)
  sim_a <- netsim(list(est_hh, est_com), param, init, control)
  sim_b <- netsim(list(est_hh, est_com), param, init, control)
  sim <- merge(sim_a, sim_b)
  expect_equal(sim$control$nsims, 2)
  expect_s3_class(sim$nwparam[[1]], "netclique")
  expect_equal(ncol(sim$epi$num), 2)
  expect_length(sim$stats$transmat, 2)
  tsim <- truncate_sim(sim, at = 5)
  expect_equal(nrow(tsim$epi$num), 6)
})
