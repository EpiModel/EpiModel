context("Network estimation functions")


# Test formation models ---------------------------------------------------

test_that("netest works for edges only model", {
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")
})


test_that("netest works for edges + nodematch model", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")
})


test_that("netest works with offset.coef terms", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "role", rep(c("I", "V", "R"), c(10, 80, 10)))
  est <- netest(nw, formation = ~edges + offset(nodematch("role", diff = TRUE,
                                                          levels = 1:2)),
                coef.form = c(-Inf, -Inf), target.stats = c(40),
                coef.diss = dissolution_coefs(~offset(edges), 52 * 2, 0.0009),
                verbose = FALSE)
  expect_is(est, "netest")
})



# Test dissolution models -------------------------------------------------

test_that("netest works for heterogeneous dissolutions", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(50, 20),
                coef.diss = dissolution_coefs(~offset(edges) + offset(nodematch("race")),
                                              c(10, 20)),
                verbose = FALSE
  )
  expect_is(est, "netest")
})

test_that("netest diss_check flags bad models", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

  formation <- ~edges + nodematch("race", diff = TRUE)
  dissolution <- ~offset(edges) + offset(nodematch("race"))
  coef.diss <- dissolution_coefs(dissolution, c(10, 20))
  expect_error(netest(nw, formation, target.stats = c(50, 20), coef.diss,
                      verbose = FALSE),
               "Term options for one or more terms in dissolution model")
})


# Other tests ---------------------------------------------------------------

test_that("Error if incorrect coef.diss parameter", {
  nw <- network_initialize(n = 50)
  coef.diss <- 1
  expect_error(netest(nw, formation = ~edges, target.stats = 25,
                      coef.diss = coef.diss,
                verbose = FALSE))
})

test_that("update_dissolution tests", {
  nw <- network_initialize(n = 1000)

  diss.300 <- dissolution_coefs(~offset(edges), 300, 0.001)
  diss.200 <- dissolution_coefs(~offset(edges), 200, 0.001)

  est300 <- netest(nw = nw,
                  formation = ~edges,
                  target.stats = 500,
                  coef.diss = diss.300)

  est200 <- netest(nw = nw,
                  formation = ~edges,
                  target.stats = 500,
                  coef.diss = diss.200)

  est200.compare <- update_dissolution(est300, diss.200)
  expect_true(round(as.numeric(est200$coef.form), 3) ==
                round(as.numeric(est200.compare$coef.form), 3))

  expect_error(update_dissolution(1, diss.200),
               "old.netest must be an object of class netest")
  expect_error(update_dissolution(est300, 1),
               "new.coef.diss must be an object of class disscoef")
  est300$edapprox <- FALSE
  expect_error(update_dissolution(est300, diss.200),
               "Edges dissolution approximation must be used")
})

test_that("differing length update_dissolution tests", {
  nw <- network_initialize(n = 1000)
  nw %v% "race" <- rep(letters[1:5], length.out = 1000)

  diss_het <- dissolution_coefs(~offset(edges) +
                                  offset(nodematch("race", diff = TRUE)),
                                c(300, 200, 100, 50, 150, 225), 0.001)
  diss_hom <- dissolution_coefs(~offset(edges), 200, 0.001)

  est_het <- netest(nw = nw,
                    formation = ~edges + nodematch("race", diff = TRUE),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_het)

  est_hom <- netest(nw = nw,
                    formation = ~edges + nodematch("race", diff = TRUE),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_hom)

  est_hom_compare <- update_dissolution(est_het, diss_hom)
  est_het_compare <- update_dissolution(est_hom, diss_het)

  expect_true(all(round(as.numeric(est_hom$coef.form), 3) ==
                round(as.numeric(est_hom_compare$coef.form), 3)))

  expect_true(all(round(as.numeric(est_het$coef.form), 3) ==
                round(as.numeric(est_het_compare$coef.form), 3)))
})

test_that("duration 1 update_dissolution tests", {
  nw <- network_initialize(n = 1000)
  nw %v% "race" <- rep(letters[1:5], length.out = 1000)

  diss_het <- dissolution_coefs(~offset(edges) +
                                  offset(nodematch("race", diff = TRUE)),
                                c(300, 200, 100, 50, 150, 225), 0.001)

  diss_hom_1 <- dissolution_coefs(~offset(edges), 1, 0)


  est_het <- netest(nw = nw,
                    formation = ~edges + nodemix("race", levels2 = 1:5),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_het, nested.edapprox = FALSE)

  est_hom_1 <- netest(nw = nw,
                    formation = ~edges + nodemix("race", levels2 = 1:5),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_hom_1)

  est_hom_1_update <- update_dissolution(est_het, diss_hom_1)
  est_het_update <- update_dissolution(est_hom_1, diss_het, nested.edapprox = FALSE)

  expect_true(all(round(as.numeric(est_hom_1$coef.form), 3) ==
                round(as.numeric(est_hom_1_update$coef.form), 3)))

  expect_true(all(round(as.numeric(est_het$coef.form), 3) ==
                round(as.numeric(est_het_update$coef.form), 3)))

  dx <- netdx(est_het, dynamic = TRUE, nsteps = 5, verbose = FALSE)
  dx <- netdx(est_het, dynamic = FALSE, nsims = 5, verbose = FALSE)

  expect_error(dx <- netdx(est_hom_1, dynamic = TRUE, nsteps = 5),
               "Running dynamic diagnostics on a cross-sectional ERGM")
  dx <- netdx(est_hom_1, dynamic = FALSE, nsims = 5, verbose = FALSE)

  dx <- netdx(est_het_update, dynamic = TRUE, nsteps = 5, verbose = FALSE)
  dx <- netdx(est_het_update, dynamic = FALSE, nsims = 5, verbose = FALSE)

  expect_error(dx <- netdx(est_hom_1_update, dynamic = TRUE, nsteps = 5),
               "Running dynamic diagnostics on a cross-sectional ERGM")
  dx <- netdx(est_hom_1_update, dynamic = FALSE, nsims = 5, verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(r.num = 0, status.vector = rep("s", 1000))
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5,
                         resimulate.network = TRUE, verbose = FALSE,
                         nwupdate.FUN = NULL)
  mod <- netsim(est_het, param, init, control)
  mod <- netsim(est_het_update, param, init, control)
  mod <- netsim(est_hom_1, param, init, control)
  mod <- netsim(est_hom_1_update, param, init, control)
})

test_that("differing length non-nested update_dissolution tests", {
  nw <- network_initialize(n = 1000)
  nw %v% "race" <- rep(letters[1:5], length.out = 1000)
  nw %v% "age" <- rep(1:3, length.out = 1000)

  diss_het <- dissolution_coefs(~offset(edges) +
                                  offset(nodematch("age", diff = TRUE)),
                                c(300, 50, 150, 225), 0.001)
  diss_hom <- dissolution_coefs(~offset(edges), 200, 0.001)

  est_het <- netest(nw = nw,
                    formation = ~edges + nodematch("race", diff = TRUE),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_het,
                    nested.edapprox = FALSE)

  est_hom <- netest(nw = nw,
                    formation = ~edges + nodematch("race", diff = TRUE),
                    target.stats = c(500, 50, 50, 90, 30, 10),
                    coef.diss = diss_hom,
                    nested.edapprox = FALSE)

  est_hom_compare <- update_dissolution(est_het, diss_hom, nested.edapprox = FALSE)
  est_het_compare <- update_dissolution(est_hom, diss_het, nested.edapprox = FALSE)

  expect_true(all(round(as.numeric(est_hom$coef.form), 3) ==
                round(as.numeric(est_hom_compare$coef.form), 3)))

  expect_true(all(round(as.numeric(est_het$coef.form), 3) ==
                round(as.numeric(est_het_compare$coef.form), 3)))

  expect_equal(est_hom_compare$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_hom_compare$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_hom_compare$formation)),
               ~offset(edges))

  expect_equal(est_het_compare$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_het_compare$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_het_compare$formation)),
               ~offset(edges) + offset(nodematch("age", diff = TRUE)))
})

test_that("non-nested EDA", {
  nw <- network_initialize(n = 1000)
  nw %v% "race" <- rep(letters[1:5], length.out = 1000)
  nw %v% "age" <- rep(1:3, length.out = 1000)

  diss_het_1 <- dissolution_coefs(~offset(edges) +
                                    offset(nodematch("race", diff = TRUE)),
                                  c(300, 200, 100, 50, 150, 225), 0.001)
  diss_het_2 <- dissolution_coefs(~offset(edges) +
                                    offset(nodematch("age", diff = TRUE)),
                                  c(55, 45, 65, 80), 0.001)

  est_nest_1 <- netest(nw = nw,
                       formation = ~edges + nodematch("race", diff = TRUE),
                       target.stats = c(500, 50, 50, 90, 30, 10),
                       coef.diss = diss_het_1)

  est_non_1 <- netest(nw = nw,
                      formation = ~edges + nodematch("race", diff = TRUE),
                      target.stats = c(500, 50, 50, 90, 30, 10),
                      coef.diss = diss_het_1,
                      nested = FALSE)

  expect_error(est_nest_2 <- netest(nw = nw,
                                    formation = ~edges + nodematch("race", diff = TRUE),
                                    target.stats = c(500, 50, 50, 90, 30, 10),
                                    coef.diss = diss_het_2),
                                    "Term options for one or more terms in dissolution model")

  est_non_2 <- netest(nw = nw,
                      formation = ~edges + nodematch("race", diff = TRUE),
                      target.stats = c(500, 50, 50, 90, 30, 10),
                      coef.diss = diss_het_2,
                      nested = FALSE)

  est_nest_1_non_2 <- update_dissolution(est_nest_1, diss_het_2, nested = FALSE)
  est_non_1_non_2 <- update_dissolution(est_non_1, diss_het_2, nested = FALSE)
  expect_error(est_nest_1_nest_2 <- update_dissolution(est_nest_1, diss_het_2, nested = TRUE),
               "Term options for one or more terms in dissolution model")
  expect_error(est_non_1_nest_2 <- update_dissolution(est_non_1, diss_het_2, nested = TRUE),
               "Term options for one or more terms in dissolution model")
  est_non_2_nest_1 <- update_dissolution(est_non_2, diss_het_1, nested = TRUE)
  est_non_2_non_1 <- update_dissolution(est_non_2, diss_het_1, nested = FALSE)


  expect_true(all(round(as.numeric(est_nest_1$coef.form), 3) ==
                round(as.numeric(est_non_2_nest_1$coef.form), 3)))

  expect_true(all(round(as.numeric(est_non_1$coef.form), 3) ==
                round(as.numeric(est_non_2_non_1$coef.form), 3)))

  expect_true(all(round(as.numeric(est_non_2$coef.form), 3) ==
                round(as.numeric(est_nest_1_non_2$coef.form), 3)))

  expect_true(all(round(as.numeric(est_non_2$coef.form), 3) ==
                round(as.numeric(est_non_1_non_2$coef.form), 3)))

  expect_equal(est_nest_1$formation, ~edges + nodematch("race", diff = TRUE))

  expect_equal(est_non_1$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_non_1$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_non_1$formation)),
               ~offset(edges) + offset(nodematch("race", diff = TRUE)))

  expect_equal(est_non_2$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_non_2$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_non_2$formation)),
               ~offset(edges) + offset(nodematch("age", diff = TRUE)))

  expect_equal(est_nest_1_non_2$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_nest_1_non_2$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_nest_1_non_2$formation)),
               ~offset(edges) + offset(nodematch("age", diff = TRUE)))

  expect_equal(est_non_1_non_2$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_non_1_non_2$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_non_1_non_2$formation)),
               ~offset(edges) + offset(nodematch("age", diff = TRUE)))

  expect_equal(est_non_2_nest_1$formation, ~edges + nodematch("race", diff = TRUE))

  expect_equal(est_non_2_non_1$formation, ~Passthrough(formation) + Passthrough(dissolution))
  expect_equal(eval(quote(formation), envir = environment(est_non_2_non_1$formation)),
               ~edges + nodematch("race", diff = TRUE))
  expect_equal(eval(quote(dissolution), envir = environment(est_non_2_non_1$formation)),
               ~offset(edges) + offset(nodematch("race", diff = TRUE)))
})

test_that("environment handling in non-nested EDA", {
  nw <- network_initialize(n = 1000)
  nw %v% "race" <- rep(letters[1:5], length.out = 1000)
  nw %v% "age" <- rep(1:3, length.out = 1000)

  make_formula <- function() {
    a <- "age"
    ff <- ~offset(edges) + offset(nodematch(a, diff = TRUE))
    ff
  }

  ff <- make_formula()

  coef_diss <- dissolution_coefs(ff, c(55, 45, 65, 80))

  r <- "race"

  netest_1 <- netest(nw = nw,
                     formation = ~edges + nodematch(r, diff = TRUE),
                     target.stats = c(500, 50, 50, 90, 40, 60),
                     coef.diss = coef_diss,
                     nested.edapprox = FALSE)

  netdx_1 <- netdx(netest_1, nsims = 2, nsteps = 5, dynamic = TRUE, verbose = FALSE)

  make_formula_2 <- function() {
    x <- "race"
    ff <- ~offset(edges) + offset(nodematch(x, diff = TRUE))
    ff
  }

  ff_2 <- make_formula_2()

  coef_diss_2 <- dissolution_coefs(ff_2, c(30, 40, 50, 25, 35, 55))

  expect_error(netest_2 <- update_dissolution(netest_1, coef_diss_2,
                                              nested.edapprox = TRUE),
               "Term options for one or more terms in dissolution model")
  netest_2 <- update_dissolution(netest_1, coef_diss_2, nested.edapprox = FALSE)
  netdx_2 <- netdx(netest_2, nsims = 2, nsteps = 5, dynamic = TRUE, verbose = FALSE)
})

test_that("non-nested EDA produces expected statistic names, with or without trimming", {
  for (trim in c(FALSE, TRUE)) {
    nw <- network.initialize(10, directed = FALSE)
    nw %v% "race" <- rep(1:3, length.out = 10)
    nw %v% "age" <- rep(1:2, length.out = 10)

    coef_diss <- dissolution_coefs(~offset(edges) + offset(nodematch("age", diff = TRUE)),
                                   c(10, 25, 20))
    est <- netest(nw,
                  formation = ~edges + nodematch("race", diff = TRUE),
                  target.stats = c(10, 3, 2, 1),
                  coef.diss = coef_diss,
                  nested.edapprox = FALSE)

    if (trim == TRUE) {
      est <- trim_netest(est)
    }

    formation_names <- c("edges", "nodematch.race.1", "nodematch.race.2",
                         "nodematch.race.3", "offset(edges)", "offset(nodematch.age.1)",
                         "offset(nodematch.age.2)")

    dxs <- netdx(est, nsims = 100, dynamic = FALSE)
    expect_equal(rownames(dxs$stats.table.formation), formation_names)

    dxd <- netdx(est, nsims = 2, nsteps = 7, dynamic = TRUE)
    expect_equal(rownames(dxd$stats.table.formation), formation_names)

    param <- param.net(inf.prob = 0.3, act.rate = 0.5)
    init <- init.net(i.num = 3)
    control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
    sim <- netsim(est, param, init, control)
    expect_equal(colnames(get_nwstats(sim, mode = "list")[[1]]), formation_names)
  }
})

test_that("trimming non-nested EDA fails when it should", {
  nw <- network.initialize(10, directed = FALSE)
  nw %v% "race" <- rep(1:3, length.out = 10)
  nw %v% "age" <- rep(1:2, length.out = 10)

  coef_diss <- dissolution_coefs(~offset(edges) + offset(nodematch("age", diff = TRUE)),
                                 c(10, 25, 20))

  r <- "race"

  est <- netest(nw,
                formation = ~edges + nodematch(r, diff = TRUE),
                target.stats = c(10, 3, 2, 1),
                coef.diss = coef_diss,
                nested.edapprox = FALSE)

  est <- trim_netest(est)

  expect_error(dxs <- netdx(est, nsims = 100, dynamic = FALSE), "object 'r' not found")
  expect_error(dxd <- netdx(est, nsims = 2, nsteps = 7, dynamic = TRUE), "object 'r' not found")

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 3)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  expect_error(sim <- netsim(est, param, init, control), "object 'r' not found")

  a <- "age"
  coef_diss <- dissolution_coefs(~offset(edges) + offset(nodematch(a, diff = TRUE)),
                                 c(10, 25, 20))

  est <- netest(nw,
                formation = ~edges + nodematch("race", diff = TRUE),
                target.stats = c(10, 3, 2, 1),
                coef.diss = coef_diss,
                nested.edapprox = FALSE)

  est <- trim_netest(est)

  expect_error(dxs <- netdx(est, nsims = 100, dynamic = FALSE), "object 'a' not found")
  expect_error(dxd <- netdx(est, nsims = 2, nsteps = 7, dynamic = TRUE), "object 'a' not found")

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 3)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  expect_error(sim <- netsim(est, param, init, control), "object 'a' not found")
})

test_that("non-nested EDA with substitutions", {
  nw <- network_initialize(n = 100)
  nw %v% "race" <- rep(letters[1:3], length.out = 100)
  nw %v% "age" <- rep(1:2, length.out = 100)

  make_form <- function() {
    a <- "race"
    ff <- ~edges + nodematch(a, diff = TRUE)
    ff
  }

  make_diss_1 <- function() {
    a <- "age"
    ff <- ~offset(edges) + offset(nodematch(a, diff = TRUE))
    ff
  }

  make_diss_2 <- function() {
    a <- "race"
    ff <- ~offset(edges) + offset(nodematch(a, diff = TRUE))
    ff
  }

  make_diss_3 <- function() {
    b <- "age"
    ff <- ~offset(edges) + offset(nodematch(b, diff = TRUE))
    ff
  }

  make_diss_4 <- function() {
    b <- "race"
    ff <- ~offset(edges) + offset(nodematch(b, diff = TRUE))
    ff
  }

  formation_names_r <- c("edges", "nodematch.race.a", "nodematch.race.b",
                         "nodematch.race.c")

  formation_names_ra <- c("edges", "nodematch.race.a", "nodematch.race.b",
                          "nodematch.race.c", "offset(edges)", "offset(nodematch.age.1)",
                          "offset(nodematch.age.2)")

  formation_names_rr <- c("edges", "nodematch.race.a", "nodematch.race.b",
                          "nodematch.race.c", "offset(edges)", "offset(nodematch.race.a)",
                          "offset(nodematch.race.b)", "offset(nodematch.race.c)")

  diss_1 <- dissolution_coefs(make_diss_1(), c(10, 5, 7))
  diss_2 <- dissolution_coefs(make_diss_2(), c(5, 4, 7, 8))
  diss_3 <- dissolution_coefs(make_diss_3(), c(20, 50, 19))
  diss_4 <- dissolution_coefs(make_diss_4(), c(3, 2, 9, 4))

  formation <- make_form()

  est_1 <- netest(nw = nw,
                  formation = formation,
                  target.stats = c(100, 10, 20, 30),
                  coef.diss = diss_1,
                  nested.edapprox = FALSE)

  est_2 <- netest(nw = nw,
                  formation = formation,
                  target.stats = c(100, 10, 20, 30),
                  coef.diss = diss_2)

  run_sims <- function(est, expected_names) {
    dxs <- netdx(est, nsims = 10, dynamic = FALSE)
    dxd <- netdx(est, nsims = 2, nsteps = 5, dynamic = TRUE)
    param <- param.net(inf.prob = 0.3, act.rate = 0.5)
    init <- init.net(i.num = 3)
    control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
    sim <- netsim(est, param, init, control)

    expect_equal(rownames(dxs$stats.table.formation), expected_names)
    expect_equal(rownames(dxd$stats.table.formation), expected_names)
    expect_equal(colnames(get_nwstats(sim, mode = "list")[[1]]), expected_names)
  }

  run_sims(est_1, formation_names_ra)
  run_sims(est_2, formation_names_r)

  est_11 <- update_dissolution(est_1, diss_1, nested.edapprox = FALSE)
  est_12 <- update_dissolution(est_1, diss_2, nested.edapprox = TRUE)
  est_13 <- update_dissolution(est_1, diss_3, nested.edapprox = FALSE)
  est_14 <- update_dissolution(est_1, diss_4, nested.edapprox = FALSE)

  run_sims(est_11, formation_names_ra)
  run_sims(est_12, formation_names_r)
  run_sims(est_13, formation_names_ra)
  run_sims(est_14, formation_names_rr)

  est_21 <- update_dissolution(est_2, diss_1, nested.edapprox = FALSE)
  est_22 <- update_dissolution(est_2, diss_2, nested.edapprox = TRUE)
  est_23 <- update_dissolution(est_2, diss_3, nested.edapprox = FALSE)
  est_24 <- update_dissolution(est_2, diss_4, nested.edapprox = FALSE)

  run_sims(est_21, formation_names_ra)
  run_sims(est_22, formation_names_r)
  run_sims(est_23, formation_names_ra)
  run_sims(est_24, formation_names_rr)
})
