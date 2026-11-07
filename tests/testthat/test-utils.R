context("Utility Functions")

test_that("brewer_ramp", {

  expect_true(length(brewer_ramp(100, plt = "Spectral")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Spectral",
                                 delete.lights = FALSE)) == 100)
  expect_false(length(brewer_ramp(100, plt = "Spectral")) == 50)
  expect_false(length(brewer_ramp(100, plt = "Spectral",
                                  delete.lights = FALSE)) == 50)

  expect_true(length(brewer_ramp(100, plt = "Accent")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Accent",
                                 delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Oranges")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Oranges",
                                 delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Set1")) == 100)

  expect_error(brewer_ramp(100, plt = "Jimmy"))
  expect_error(brewer_ramp(100, plt = "Jimmy", delete.lights = FALSE))
  expect_error(brewer_ramp(-1, plt = "Spectral"))

})


test_that("color_tea", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, rec.rate = 0.001)
  init <- init.net(i.num = 25, r.num = 25)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  nd <- get_network(mod)
  nd <- color_tea(nd, verbose = FALSE)
  expect_is(nd, "networkDynamic")
  expect_true(length(unique(
    get.vertex.attribute.active(nd, "ndtvcol", at = 1))) == 3)
})

test_that("delete_attr", {
  dat <- create_dat_object()
  dat <- append_core_attr(dat, 1, 5)
  dat <- append_attr(dat, "a", 1:5, 5)
  dat <- append_attr(dat, "b", 6:10, 5)

  dat <- delete_attr(dat, 5)
  expect_is(dat, "list")

  dat <- delete_attr(dat, 2:3)
  expect_true(length(unique(vapply(get_attr_list(dat), length, 1))) == 1)

  l2 <- dat
  l3 <- delete_attr(dat, NULL)

  expect_equal(l2, l3)
})

test_that("ssample", {

  expect_true(length(ssample(1:5, 1)) == 1)
  expect_equal(ssample(5, 1), 5)
  expect_equal(ssample(5, 2), 5)
  expect_null(ssample(5, 0))

})

test_that("apportion_lr happy path returns a vector of the requested length", {
  out <- apportion_lr(20, c("a", "b", "c", "d"), c(0.5, 0.25, 0.125, 0.125))
  expect_length(out, 20)
  expect_setequal(unique(out), c("a", "b", "c", "d"))

  # Length-(k-1) proportions: missing entry fills to sum to 1.
  out2 <- apportion_lr(8, c("a", "b", "c", "d"), c(0.25, 0.25, 0.25))
  expect_length(out2, 8)

  out3 <- apportion_lr(10, 1:4, c(0.25, 0.25, 0.25, 0.25), shuffled = TRUE)
  expect_length(out3, 10)
})

test_that("apportion_lr validates its arguments", {
  expect_error(apportion_lr(1.5, 1:3, c(0.3, 0.3, 0.4)),
               "vector.length must be a positive integer")
  expect_error(apportion_lr(0, 1:3, c(0.3, 0.3, 0.4)),
               "vector.length must be a positive integer")
  # is.vector() returns TRUE for lists, so we need a non-vector type
  # (matrix here) to trip the "values must be a vector" guard.
  expect_error(apportion_lr(5, matrix(1:6, 2), c(0.3, 0.3, 0.4)),
               "values must be a vector")
  # Length matches but sum != 1 and not within the (length - 1, [0, 1]) branch.
  expect_error(apportion_lr(5, 1:3, c(0.3, 0.3, 0.3)),
               "proportions length or proportions sum")
})

test_that("check_degdist_bal", {
  expect_output(check_degdist_bal(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g2 = c(0.40, 0.55, 0.03, 0.02),
                                  deg.dist.g1 = c(0.48, 0.41, 0.08, 0.03)),
                "-0.015 Rel Diff")
  expect_output(check_degdist_bal(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03)),
                "Edges balanced")
  expect_output(check_degdist_bal(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.45, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03)),
                "deg.dist.g1 TOTAL != 1")
  expect_output(check_degdist_bal(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.55, 0.41, 0.08, 0.03)),
                "deg.dist.g2 TOTAL != 1")
})

test_that("cap_ncores handles unknown detected core counts", {
  expect_equal(EpiModel:::cap_ncores(2, detected = NA_integer_), 2)
  expect_equal(EpiModel:::cap_ncores(2, detected = 0), 2)
  expect_equal(EpiModel:::cap_ncores(2, detected = Inf), 2)
  expect_equal(EpiModel:::cap_ncores(2, detected = 4), 2)
  expect_equal(EpiModel:::cap_ncores(4, detected = 2), 2)
  expect_error(EpiModel:::cap_ncores(NA_real_, detected = 2), "`ncores` must be")
})

test_that("edgelist_censor", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  sim <- netdx(est, nsims = 1, nsteps = 100,
               keep.tedgelist = TRUE, verbose = FALSE)
  el <- as.data.frame(sim)
  expect_is(edgelist_censor(el), "matrix")
})

test_that("get_degree", {
  nw <- network_initialize(n = 500)

  set.seed(1)
  fit <- ergm(nw ~ edges, target.stats = 250, eval.loglik = FALSE)
  sim <- simulate(fit)

  # Slow ERGM-based method
  ergm.method <- unname(summary(sim ~ sociality(nodes = TRUE)))
  ergm.method

  # Fast tabulate method with network object
  deg.net <- get_degree(sim)
  deg.net

  # Even faster if network already transformed into an edgelist
  el <- as.edgelist(sim)
  deg.el <- get_degree(el)
  deg.el

  expect_equal(ergm.method, deg.net)
  expect_equal(ergm.method, deg.el)
})

test_that("dissolution_coefs stops at the departure ceiling and names the duration", {
  # Ceiling is 1 / (1 - (1 - d.rate)^2): 30.25 at d.rate = 1/60, 50.25 at 0.01
  expect_error(dissolution_coefs(~offset(edges), duration = 60, d.rate = 1 / 60),
               "^duration = 60 is at or above the departure ceiling of 30.2521 for d.rate = 0.01666667")
  expect_error(dissolution_coefs(~offset(edges), duration = 60, d.rate = 0.01),
               "ceiling of 50.2513 for d.rate = 0.01 .*Lower the duration below the ceiling, or if it is right, lower d.rate below 0.008368\\.$")
  # The stated d.rate bound admits the duration, which then sits at the ceiling
  expect_warning(o <- dissolution_coefs(~offset(edges), duration = 60, d.rate = 0.008368),
                 "duration = 60 \\(100.0%\\) is at or above 90% of the departure ceiling of 60.0025")
  expect_s3_class(o, "disscoef")
  expect_error(dissolution_coefs(~offset(edges), duration = 60, d.rate = 0.008369),
               "at or above the departure ceiling")

  dissolution <- ~offset(edges) + offset(nodematch("age.grp", diff = TRUE))
  duration <- c(60, 30, 80, 100, 125, 160)
  expect_error(dissolution_coefs(dissolution, duration = duration, d.rate = 1 / 60),
               "^duration\\[1\\] = 60, duration\\[3\\] = 80, duration\\[4\\] = 100, duration\\[5\\] = 125 and duration\\[6\\] = 160 are at or above the departure ceiling of 30.2521 .*Lower these durations below the ceiling, or if they are right, lower d.rate below 0.003129\\.$")
  # The error names the dissolution_coefs call, not the internal check
  err <- tryCatch(dissolution_coefs(~offset(edges), duration = 60, d.rate = 0.01),
                  error = function(e) e)
  expect_equal(conditionCall(err)[[1]], as.name("dissolution_coefs"))
})

test_that("dissolution_coefs warns at or above 90% of the departure ceiling", {
  d.rate <- 0.000478213
  ceiling <- 1 / (1 - (1 - d.rate)^2)
  expect_equal(ceiling, 1045.809, tolerance = 1e-6)

  expect_warning(o <- dissolution_coefs(~offset(edges), duration = 1040, d.rate = d.rate),
                 "^duration = 1040 \\(99.4%\\) is at or above 90% of the departure ceiling of 1045.81 for d.rate = 0.000478213 .*Check the duration first; d.rate is the other lever\\.$")
  # The coefficients are still returned, and the adjustment is -log(1 - share)
  expect_equal(o$coef.crude, log(1039))
  expect_equal(o$coef.adj, o$coef.crude - log(1 - 1040 / ceiling))
  expect_equal(o$coef.adj, 12.139, tolerance = 1e-4)
  warn <- tryCatch(dissolution_coefs(~offset(edges), duration = 1040, d.rate = d.rate),
                   warning = function(w) w)
  expect_equal(conditionCall(warn)[[1]], as.name("dissolution_coefs"))

  # Heterogeneous model: one warning naming every duration at or above 90%
  dissolution <- ~offset(edges) + offset(nodematch("age.grp", diff = TRUE))
  expect_warning(o <- dissolution_coefs(dissolution, duration = c(60, 30, 950, 100, 1040, 160),
                                        d.rate = d.rate),
                 "^duration\\[3\\] = 950 \\(90.8%\\) and duration\\[5\\] = 1040 \\(99.4%\\) are at or above 90% of the departure ceiling of 1045.81 .*Check these durations first")
  expect_length(o$coef.adj, 6)

  # Below 90% of the ceiling, and with no departures, there is no warning
  expect_no_warning(dissolution_coefs(~offset(edges), duration = 900, d.rate = d.rate))
  expect_no_warning(dissolution_coefs(~offset(edges), duration = 1e6, d.rate = 0))
  expect_no_warning(dissolution_coefs(dissolution, duration = c(60, 30, 900, 100, 125, 160),
                                      d.rate = d.rate))
})

test_that("dissolution_coefs coefficients are unchanged", {
  o <- dissolution_coefs(~offset(edges), duration = 38, d.rate = 0.01)
  expect_equal(o$coef.crude, 3.61091791264423)
  expect_equal(o$coef.adj, 5.0223249745792)
  expect_equal(o$coef.form.corr, 3.63758615972639)
  expect_equal(o$diss.model.type, "edgesonly")

  o <- dissolution_coefs(~offset(edges) + offset(nodematch("race")), duration = c(20, 10),
                         d.rate = 0.001)
  expect_equal(o$coef.crude, c(2.94443897916644, -0.74721440183022))
  expect_equal(o$coef.adj, c(2.98524014057037, -0.767823059946204))
  expect_equal(o$coef.form.corr, c(2.99573227355399, -0.693147180559944))
  expect_equal(o$diss.model.type, "nodematch")

  o <- dissolution_coefs(~offset(edges) + offset(nodemix("race")),
                         duration = c(60, 30, 80, 100, 125, 160), d.rate = 0.002)
  expect_equal(o$coef.crude, c(4.07753744390572, -0.710241613919242, 0.291910408561309,
                               0.517582406228873, 0.74274412169932, 0.991366758314519))
  expect_equal(o$coef.adj, c(4.35165854998479, -0.856665702828121, 0.402981305750843,
                             0.753620479372626, 1.1607706958471, 1.73712070036617))
  expect_equal(o$coef.form.corr, c(4.0943445622221, -0.693147180559942, 0.287682072451787,
                                   0.510825623765993, 0.733969175080203, 0.980829253011732))
  expect_equal(o$diss.model.type, "nodemix")

  expect_error(dissolution_coefs(~offset(edges), duration = c(20, 10)),
               "Dissolution model length is 1, but number of duration was 2")
  expect_error(dissolution_coefs(~offset(edges) + offset(nodematch("a")) + offset(nodemix("b")),
                                 duration = c(20, 10, 5)),
               "The form of the dissolution argument is invalid")
})


test_that("get_formula_term_attr checks", {
  nw <- network_initialize(n = 100)

  expect_null(get_formula_term_attr(~edges, nw))

  nw <- network_initialize(n = 100)
  riskg <- sample(rep(1:2, each = 50))
  race <- sample(rep(0:1, each = 50))
  nw <- set_vertex_attribute(nw, "riskg", riskg)
  nw <- set_vertex_attribute(nw, "race", race)
  expect_null(get_formula_term_attr(~edges, nw))

  expect_equal(get_formula_term_attr(~edges + nodefactor("race"), nw), "race")
  expect_equal(get_formula_term_attr(~edges + nodefactor("race") +
                                       nodematch("riskg"), nw),
               c("race", "riskg"))
  expect_equal(get_formula_term_attr(~edges + nodefactor(c("race", "riskg")) +
                                       nodematch("riskg"), nw),
               c("race", "riskg"))

})

context("Shiny App")

test_that("epiweb works", {
  expect_error(epiweb(class = "foo"))
})

context("update.R Functionality")

test_that("add_vertices, delete_vertices, delete_edges behave as expected", {
  net_size <- 10L
  nw <- network_initialize(net_size)
  nw[1, 2] <- 1
  nw[3, 6] <- 1
  nw[2, 5] <- 1
  nw[8, 10] <- 1
  el <- as.edgelist(nw)
  el <- matrix(c(el), ncol = 2)
  attr(el, "n") <- net_size

  ## add some vertices
  el_add <- add_vertices(el, 4)
  expect_equal(el, el_add, check.attributes = FALSE)
  expect_equal(attr(el_add, "n"), net_size + 4)

  ## delete no vertices, no edges
  vertices <- NULL
  el_del <- matrix(c(1, 2, 3, 8, 2, 5, 6, 10), ncol = 2)
  el_drop <- matrix(c(1, 2, 3, 8, 2, 5, 6, 10), ncol = 2)
  el_del_edges <- delete_edges(el, vertices)
  el_del_verts <- delete_vertices(el, vertices)
  expect_equal(el_del_edges, el_del, check.attributes = FALSE)
  expect_equal(el_del_verts, el_drop, check.attributes = FALSE)
  expect_equal(attr(el_del_edges, "n"), net_size)
  expect_equal(attr(el_del_verts, "n"), net_size - length(vertices))

  ## delete some vertices, no edges
  vertices <- c(4, 7)
  el_del <- matrix(c(1, 2, 3, 8, 2, 5, 6, 10), ncol = 2)
  el_drop <- matrix(c(1, 2, 3, 6, 2, 4, 5, 8), ncol = 2)
  el_del_edges <- delete_edges(el, vertices)
  el_del_verts <- delete_vertices(el, vertices)
  expect_equal(el_del_edges, el_del, check.attributes = FALSE)
  expect_equal(el_del_verts, el_drop, check.attributes = FALSE)
  expect_equal(attr(el_del_edges, "n"), net_size)
  expect_equal(attr(el_del_verts, "n"), net_size - length(vertices))

  ## delete some vertices, some edges
  vertices <- c(2, 8, 9)
  el_del <- matrix(c(3, 6), ncol = 2)
  el_drop <- matrix(c(2, 5), ncol = 2)
  el_del_edges <- delete_edges(el, vertices)
  el_del_verts <- delete_vertices(el, vertices)
  expect_equal(el_del_edges, el_del, check.attributes = FALSE)
  expect_equal(el_del_verts, el_drop, check.attributes = FALSE)
  expect_equal(attr(el_del_edges, "n"), net_size)
  expect_equal(attr(el_del_verts, "n"), net_size - length(vertices))

  ## delete some vertices, all edges
  vertices <- c(1, 2, 3, 4, 5, 6, 8, 10)
  el_del <- matrix(integer(0), ncol = 2)
  el_drop <- matrix(integer(0), ncol = 2)
  el_del_edges <- delete_edges(el, vertices)
  el_del_verts <- delete_vertices(el, vertices)
  expect_equal(el_del_edges, el_del, check.attributes = FALSE)
  expect_equal(el_del_verts, el_drop, check.attributes = FALSE)
  expect_equal(attr(el_del_edges, "n"), net_size)
  expect_equal(attr(el_del_verts, "n"), net_size - length(vertices))

  ## delete all vertices, all edges
  vertices <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  el_del <- matrix(integer(0), ncol = 2)
  el_drop <- matrix(integer(0), ncol = 2)
  el_del_edges <- delete_edges(el, vertices)
  el_del_verts <- delete_vertices(el, vertices)
  expect_equal(el_del_edges, el_del, check.attributes = FALSE)
  expect_equal(el_del_verts, el_drop, check.attributes = FALSE)
  expect_equal(attr(el_del_edges, "n"), net_size)
  expect_equal(attr(el_del_verts, "n"), net_size - length(vertices))
})

test_that("shiftVec matches its R reference and handles edge cases", {
  # x minus the number of deleted ids at or below it
  shift_ref <- function(x, y) as.integer(x - findInterval(x, y))
  set.seed(1)
  for (k in 1:50) {
    x <- sort(sample.int(1000, sample(0:500, 1), replace = TRUE))
    y <- sort(sample.int(1000, sample(0:100, 1), replace = TRUE))
    expect_identical(EpiModel:::shiftVec(x, y), shift_ref(x, y))
  }
  # double input is coerced and the result is integer
  expect_identical(EpiModel:::shiftVec(c(1, 5, 9), c(3, 6)), c(1L, 4L, 7L))
  expect_identical(EpiModel:::shiftVec(integer(0), 1:3), integer(0))
  expect_identical(EpiModel:::shiftVec(1:5, integer(0)), 1:5)
  # deleted ids all below, all above, and equal to ids in x
  expect_identical(EpiModel:::shiftVec(c(5L, 6L), c(1L, 2L)), c(3L, 4L))
  expect_identical(EpiModel:::shiftVec(c(5L, 6L), c(8L, 9L)), c(5L, 6L))
  expect_identical(EpiModel:::shiftVec(c(3L, 5L), c(3L, 3L)), c(1L, 3L))
  expect_error(EpiModel:::shiftVec("a", 1L), "x and y must be numeric vectors")
})
