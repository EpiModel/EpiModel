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
                                                          keep = 1:2)),
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
                coef.diss = dissolution_coefs(~offset(edges) +
                                                offset(nodematch("race")),
                                              c(10, 20)),
                verbose = FALSE
  )
  expect_is(est, "netest")
})

test_that("netest diss_check flags bad models", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

  formation <- ~edges + nodematch("race")
  dissolution <- ~offset(edges) + offset(nodefactor("race"))
  coef.diss <- dissolution_coefs(dissolution, c(10, 20))
  expect_error(netest(nw, formation, target.stats = c(50, 20),
               coef.diss, verbose = FALSE),
               "Dissolution model is not a subset of formation model.")

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

  expect_equal(est_hom_compare$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges))
  expect_equal(est_het_compare$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("age", diff = TRUE)))
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
  expect_equal(est_non_1$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("race", diff = TRUE)))
  expect_equal(est_non_2$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("age", diff = TRUE)))

  expect_equal(est_nest_1_non_2$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("age", diff = TRUE)))
  expect_equal(est_non_1_non_2$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("age", diff = TRUE)))
  expect_equal(est_non_2_nest_1$formation, ~edges + nodematch("race", diff = TRUE))
  expect_equal(est_non_2_non_1$formation, ~edges + nodematch("race", diff = TRUE) +
                 offset(edges) + offset(nodematch("race", diff = TRUE)))
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

  expect_error(netdx_1 <- netdx(netest_1, nsims = 2, nsteps = 5, dynamic = TRUE),
               "object 'a' not found")

  netest_2 <- netest(nw = nw,
                     formation = ~edges + nodematch(r, diff = TRUE),
                     target.stats = c(500, 50, 50, 90, 40, 60),
                     coef.diss = coef_diss,
                     nested.edapprox = FALSE,
                     from.new = "a")

  expect_error(netdx_2 <- netdx(netest_2, nsims = 2, nsteps = 5, dynamic = TRUE), NA)

  make_formula_2 <- function() {
    x <- "race"
    ff <- ~offset(edges) + offset(nodematch(x, diff = TRUE))
    ff
  }

  ff_2 <- make_formula_2()

  coef_diss_2 <- dissolution_coefs(ff_2, c(30, 40, 50, 25, 35, 55))

  expect_error(netest_3 <- update_dissolution(netest_2, coef_diss_2,
                                              nested.edapprox = TRUE),
               "Term options for one or more terms in dissolution model")
  netest_3 <- update_dissolution(netest_2, coef_diss_2, nested.edapprox = FALSE)
  expect_error(netdx_3 <- netdx(netest_3, nsims = 2, nsteps = 5, dynamic = TRUE),
               "object 'x' not found")
  netest_4 <- update_dissolution(netest_2, coef_diss_2, nested.edapprox = FALSE,
                                 from.new = "x")
  expect_error(netdx_4 <- netdx(netest_4, nsims = 2, nsteps = 5, dynamic = TRUE), NA)
})

# STERGM --------------------------------------------------------------------

test_that("Basic STERGM fit", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                edapprox = FALSE, verbose = FALSE)
  expect_is(est, "netest")
  expect_true(!est$edapprox)
})
