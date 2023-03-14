context("transmat functions")

# test that transmat output has right class
require(EpiModel)
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 100, nsims = 1, verbose = FALSE)
mod1 <- netsim(est1, param, init, control)
tm <- get_transmat(mod1)

test_that("transmat class", {
  expect_true(is.transmat(tm))
  expect_false(is.transmat(list()))
  expect_true(inherits(tm, c("transmat", "data.frame")))
  expect_equal(names(tm), c("at", "sus", "inf", "network", "infDur", "transProb",
                            "actRate", "finalProb"))
})

# define a known transmat output
tm <- structure(list(at = c(2L, 2L, 4L, 6L, 10L, 11L, 11L, 12L, 12L,  13L, 13L,
                            13L, 13L, 14L, 14L, 15L, 15L, 16L, 17L, 19L, 19L,
                            21L,  22L, 22L, 23L, 24L, 25L, 25L, 25L, 26L, 26L,
                            26L, 26L, 27L, 27L,  31L, 31L, 31L, 32L, 32L, 32L,
                            32L, 35L, 35L, 36L, 36L, 37L, 37L,  37L, 39L, 40L,
                            40L, 40L),
                     sus = c(83L, 24L, 75L, 5L, 76L, 9L,  8L, 10L, 38L, 36L,
                             35L, 21L, 97L, 87L, 66L, 39L, 17L, 78L, 98L, 92L,
                             82L, 70L, 81L, 60L, 15L, 18L, 51L, 90L, 28L, 14L,
                             45L, 99L,  56L, 7L, 73L, 27L, 49L, 64L, 46L, 62L,
                             68L, 42L, 44L, 79L, 69L,  63L, 86L, 25L, 1L, 48L,
                             32L, 3L, 37L),
                     inf = c(72L, 72L, 72L,  83L, 83L, 83L, 76L, 8L, 8L, 10L,
                             10L, 10L, 75L, 8L, 36L, 10L,  36L, 66L, 17L, 78L,
                             36L, 92L, 83L, 70L, 60L, 24L, 39L, 70L, 18L, 70L,
                             51L, 81L, 51L, 18L, 36L, 36L, 81L, 8L, 49L, 21L,
                             14L, 8L,  62L, 17L, 38L, 44L, 69L, 73L, 69L, 99L,
                             35L, 18L, 1L),
                     infDur = c(29,  29, 31, 4, 8, 9, 1, 1, 1, 1, 1, 1, 9, 3, 1,
                                3, 2, 2, 2, 3, 6,  2, 20, 1, 1, 22, 10, 4, 1, 5,
                                1, 4, 1, 3, 14, 18, 9, 20, 1, 19,  6, 21, 3, 20,
                                24, 1, 1, 10, 1, 13, 27, 16, 3),
                     transProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                     actRate = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1),
                     finalProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)),
                .Names = c("at", "sus", "inf", "infDur", "transProb", "actRate",
                           "finalProb"),
                row.names = c(NA, 53L), class = c("transmat", "data.frame"))


test_that("as.network.transmat", {
  skip_on_cran()
  net <- as.network(tm)
  expect_equal(network.size(net), 54)
  expect_equal(list.edge.attributes(net),
               c("actRate", "at", "finalProb", "infDur", "na", "transProb"))
  expect_equal(list.vertex.attributes(net), c("at", "na", "vertex.names"))

  plot(tm) #phylogram default
  plot(tm, style = "network")
})
