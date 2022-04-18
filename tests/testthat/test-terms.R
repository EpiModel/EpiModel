context("EpiModel terms")

nw <- network_initialize(n = 50)
age <- runif(50)
sex <- rep(c(0,1), length.out = 50)
nw %v% "age" <- age
nw %v% "sex" <- sex

test_that("netest works for EpiModel terms", {
  est <- netest(nw, formation = ~edges + absdiffby("age", "sex", 1) + nodematch("sex"), target.stats = c(25, 25, 0),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")

  est2 <- netest(nw, formation = ~edges + absdiffnodemix("age", "sex"), target.stats = c(40, 5, 10, 5),
                 coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                 verbose = FALSE)
  expect_is(est2, "netest")
})

test_that("EpiModel terms produce correct summary statistics", {
  nw1 <- san(nw ~ edges + offset(nodematch("sex")), target.stats = c(30), offset.coef = c(-Inf))
    
  el1 <- as.edgelist(nw1)
  
  expect_equal(summary(nw1 ~ absdiffby("age", "sex", 2.3)),
               sum(abs(age[el1[,1]] - 2.3*sex[el1[,1]] - age[el1[,2]] + 2.3*sex[el1[,2]])),
               check.attributes = FALSE)
  
  nw2 <- san(nw ~ edges, target.stats = c(40))

  el2 <- as.edgelist(nw2)

  expect_equal(summary(nw2 ~ absdiffnodemix("age", "sex")),
               c(sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 0 & sex[el2[,2]] == 0)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(((sex[el2[,1]] == 0) + (sex[el2[,2]] == 0)) == 1)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 1 & sex[el2[,2]] == 1))),
               check.attributes = FALSE)
})

test_that("EpiModel terms produce correct change statistics", {
  nw1 <- simulate(nw ~ edges + offset(nodematch("sex")),
                  coef = c(-3, -Inf),
                  monitor = ~absdiffby("age", "sex", 2.3))
  
  stats1 <- attr(nw1, "stats")
  
  el1 <- as.edgelist(nw1)
  
  expect_equal(stats1[-c(1,2)],
               sum(abs(age[el1[,1]] - 2.3*sex[el1[,1]] - age[el1[,2]] + 2.3*sex[el1[,2]])),
               check.attributes = FALSE)
  
  nw2 <- simulate(nw ~ edges,
                  coef = c(-3),
                  monitor = ~absdiffnodemix("age", "sex"))
  
  stats2 <- attr(nw2, "stats")
  
  el2 <- as.edgelist(nw2)

  expect_equal(stats2[-c(1)],
               c(sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 0 & sex[el2[,2]] == 0)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(((sex[el2[,1]] == 0) + (sex[el2[,2]] == 0)) == 1)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 1 & sex[el2[,2]] == 1))),
               check.attributes = FALSE)
})
