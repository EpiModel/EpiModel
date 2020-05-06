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
  expect_equal(class(tm), c("transmat", "data.frame"))
  expect_equal(names(tm), c( "at", "sus", "inf", "infDur", "transProb", "actRate", "finalProb"))
})

# define a known transmat output
tm <- structure(list(at = c(2L, 2L, 4L, 6L, 10L, 11L, 11L, 12L, 12L,  13L, 13L, 13L, 13L,
                            14L, 14L, 15L, 15L, 16L, 17L, 19L, 19L, 21L,  22L, 22L, 23L,
                            24L, 25L, 25L, 25L, 26L, 26L, 26L, 26L, 27L, 27L,  31L, 31L,
                            31L, 32L, 32L, 32L, 32L, 35L, 35L, 36L, 36L, 37L, 37L,  37L,
                            39L, 40L, 40L, 40L),
                     sus = c(83L, 24L, 75L, 5L, 76L, 9L,  8L,
                            10L, 38L, 36L, 35L, 21L, 97L, 87L, 66L, 39L, 17L, 78L,
                            98L,  92L, 82L, 70L, 81L, 60L, 15L, 18L, 51L, 90L, 28L, 14L,
                            45L, 99L,  56L, 7L, 73L, 27L, 49L, 64L, 46L, 62L, 68L, 42L,
                            44L, 79L, 69L,  63L, 86L, 25L, 1L, 48L, 32L, 3L, 37L),
                     inf = c(72L, 72L, 72L,  83L, 83L, 83L, 76L, 8L, 8L, 10L, 10L, 10L,
                             75L, 8L, 36L, 10L,  36L, 66L, 17L, 78L, 36L, 92L, 83L, 70L,
                             60L, 24L, 39L, 70L, 18L,  70L, 51L, 81L, 51L, 18L, 36L, 36L,
                             81L, 8L, 49L, 21L, 14L, 8L,  62L, 17L, 38L, 44L, 69L, 73L, 69L,
                             99L, 35L, 18L, 1L),
                     infDur = c(29,  29, 31, 4, 8, 9, 1, 1, 1, 1, 1, 1, 9, 3, 1, 3, 2, 2, 2,
                                3, 6,  2, 20, 1, 1, 22, 10, 4, 1, 5, 1, 4, 1, 3, 14, 18, 9,
                                20, 1, 19,  6, 21, 3, 20, 24, 1, 1, 10, 1, 13, 27, 16, 3),
                     transProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5 ),
                     actRate = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                     finalProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                   0.5, 0.5, 0.5, 0.5, 0.5 )),
                .Names = c("at", "sus", "inf", "infDur", "transProb", "actRate",  "finalProb"),
                row.names = c(NA, 53L), class = c("transmat", "data.frame" ))

# check that it reconstructs tm into the expected phylo object
test_that("as.phylo.transmat",{
  skip_on_cran()
  phy <- as.phylo.transmat(tm)
  expect_equal(class(phy), "phylo")
  expect_equal(names(phy),
               c("edge", "Nnode", "tip.label", "node.label", "root.edge", "edge.length"))
  # lazy test, this is just veryify that the output matches previous output.
  expect_equal(phy,structure(list(edge = structure(c(55, 55, 56, 56, 57, 57, 58,  58, 59, 59, 60,
                                                     60, 61, 61, 62, 62, 63, 63, 64, 64, 65, 65, 66,
                                                     66, 67, 67, 68, 68, 69, 69, 70, 70, 71, 71, 72, 72,
                                                     73, 73, 74,  74, 75, 75, 76, 76, 77, 77, 78, 78,
                                                     79, 79, 80, 80, 81, 81, 82,  82, 83, 83, 84, 84,
                                                     85, 85, 86, 86, 87, 87, 88, 88, 89, 89, 90,  90,
                                                     91, 91, 92, 92, 93, 93, 94, 94, 95, 95, 96, 96,
                                                     97, 97, 98,  98, 99, 99, 100, 100, 101, 101, 102,
                                                     102, 103, 103, 104, 104,  105, 105, 106, 106, 107,
                                                     107, 56, 58, 57, 80, 1, 67, 59, 30,  60, 61, 77, 31,
                                                     3, 62, 63, 64, 68, 99, 65, 69, 66, 105, 70, 94,  6,
                                                     32, 92, 33, 71, 72, 5, 81, 75, 73, 8, 74, 98, 34, 10,
                                                     76,  89, 35, 11, 78, 2, 86, 82, 79, 13, 36, 14, 83,
                                                     15, 85, 84, 37,  88, 38, 12, 95, 87, 39, 91, 104, 17,
                                                     40, 106, 41, 90, 102, 7,  42, 18, 93, 96, 43, 19, 44,
                                                     20, 97, 21, 45, 4, 46, 22, 100, 9,  47, 23, 101, 24,
                                                     48, 103, 49, 26, 50, 25, 107, 27, 51, 28, 52,  16, 53,
                                                     29, 54), .Dim = c(106L, 2L),
                                                   .Dimnames = list(NULL, c("",  "infector"))), Nnode = 53L,
                                  tip.label = c(72L, 83L, 76L, 8L,  10L, 75L, 36L, 66L, 17L, 78L, 92L, 70L,
                                                60L, 24L, 39L, 18L, 51L,  81L, 49L, 21L, 14L, 62L, 38L, 44L,
                                                69L, 73L, 99L, 35L, 1L, 5L,  9L, 97L, 87L, 98L, 82L, 15L, 90L,
                                                28L, 45L, 56L, 7L, 27L, 64L,  46L, 68L, 42L, 79L, 63L, 86L, 25L,
                                                48L, 32L, 3L, 37L),
                                  node.label = c(72L,  72L, 72L, 83L, 83L, 83L, 76L, 8L, 8L, 10L, 10L, 10L, 75L,
                                                 8L,  36L, 10L, 36L, 66L, 17L, 78L, 36L, 92L, 83L, 70L, 60L, 24L,
                                                 39L,  70L, 18L, 70L, 51L, 81L, 51L, 18L, 36L, 36L, 81L, 8L, 49L,
                                                 21L,  14L, 8L, 62L, 17L, 38L, 44L, 69L, 73L, 69L, 99L, 35L, 18L, 1L),
                                  root.edge = 2L, edge.length = c(0, 4, 2, 22, 37, 9, 4, 35,  1, 1, 11, 30, 30,
                                                                  1, 0, 1, 2, 24, 0, 1, 0, 27, 2, 19, 28, 28, 17,
                                                                  27, 1, 2, 26, 10, 4, 2, 25, 3, 18, 24, 22, 2, 8,
                                                                  22, 20,  1, 19, 4, 3, 1, 18, 18, 17, 1, 16, 1, 1,
                                                                  16, 2, 16, 15, 6, 0,  15, 5, 13, 15, 15, 13, 14,
                                                                  4, 10, 10, 10, 10, 1, 1, 10, 9, 9,  9, 3, 9, 9, 9,
                                                                  9, 6, 1, 6, 6, 5, 1, 5, 5, 0, 4, 4, 4, 4, 3, 2,  2,
                                                                  1, 1, 1, 1, 1, 1)),
                             .Names = c("edge", "Nnode", "tip.label",  "node.label", "root.edge", "edge.length"),
                             class = "phylo"))

  # test method dispatch, requires ape loaded
  require(ape)
  tiny <- matrix(c(1,2,2,2,3,4,0,5,10), ncol = 3)
  tiny <- as.data.frame(tiny)
  class(tiny) <- c('transmat',class(tiny))

  # expect warning because colnames missing
  expect_warning(as.phylo(tiny), regexp = "input does not have appropriate column names")
  colnames(tiny) <- c('inf', 'sus', 'at')
  as.phylo(tiny)

  # check adding departure times
  expect_equal(as.phylo(tiny, vertex.exit.times = c(3,11,7,15))$edge.length,
               c(3,5,5,2,1,5))

  # check departure time with NA
  expect_equal(as.phylo(tiny,vertex.exit.times = c(3,11,7,NA))$edge.length,
               c(3,5,5,2,1,2))

  # more complete vertex.exit time example

  evo_tm <- structure(list(at = c(50L, 230L, 573L, 622L, 650L, 652L, 678L,  852L,
                                  914L, 1032L, 1083L, 1120L, 1134L, 1141L, 1154L, 1252L,
                                  1296L, 1320L, 1329L, 1349L, 1351L, 1354L, 1385L, 1427L,
                                  1451L,  1457L, 1474L, 1518L, 1520L, 1521L, 1533L, 1534L,
                                  1540L, 1541L,  1551L, 1566L, 1572L, 1581L, 1588L, 1603L,
                                  1622L, 1624L, 1628L,  1650L, 1651L, 1653L, 1653L, 1670L,
                                  1687L, 1711L, 1714L, 1723L,  1759L, 1800L, 1818L, 1821L),
                           sus = c(89, 66, 25, 22, 38, 58,  91, 18, 81, 40, 62, 48, 111, 76, 80,
                                   56, 53, 39, 94, 9, 102,  11, 17, 33, 6, 12, 64, 10, 3, 82,
                                   107, 108, 68, 41, 26, 78, 109,  4, 114, 106, 42, 92, 103, 30,
                                   75, 112, 93, 7, 43, 24, 95, 60,  110, 2, 55, 54),
                           inf = c(83, 83, 83, 25, 22, 22, 25, 83, 22,  22, 40, 83, 62, 40, 111,
                                   48, 66, 81, 53, 94, 62, 40, 80, 25,  33, 48, 56, 18, 81, 81,
                                   25, 56, 10, 66, 39, 107, 10, 48, 48,  9, 62, 108, 40, 11, 106,
                                   94, 103, 108, 75, 25, 18, 9, 10, 82,  83, 110)),
                      .Names = c("at", "sus", "inf"), row.names = c(NA,  -56L), class = "transmat")
  evo_exit_times <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 772,  NA, NA, NA,
                      NA, NA, NA, 1169, NA, NA, NA, 1599, NA, 223, NA,  NA, 165, NA, NA, NA,
                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA,
                      NA, NA, NA, NA, NA, 1713, 1322, NA,  NA, NA, NA, NA, NA, NA, NA, NA,
                      NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, 1090,
                      NA, NA, NA, NA, 1651, NA,  NA, NA, NA, 1796, NA, NA, NA, NA, NA, NA,
                      NA, NA, 1727, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                      NA, NA, NA, NA, NA,  NA, NA)

  evoPhylo <- as.phylo(evo_tm,vertex.exit.times = evo_exit_times)
  expect_equal(evoPhylo$edge.length, c(180, 1601, 343, 1066, 279, 49, 56, 28, 2, 1173, 262,
                                       670, 749,  1145, 268, 666, 118, 406, 137, 51, 58, 51,
                                       698, 132, 217, 20,  213, 682, 669, 231, 205, 222, 245,
                                       33, 200, 231, 494, 20, 304,  254, 271, 472, 274, 296,
                                       438, 438, 106, 24, 372, 372, 124, 366,  60, 349, 196,
                                       22, 1, 303, 302, 279, 178, 33, 289, 90, 32, 283,  282,
                                       282, 272, 48, 257, 257, 187, 251, 7, 242, 235, 235, 120,
                                       48, 201, 201, 46, 199, 195, 25, 173, 173, 172, 36, 143, 170,
                                       74, 170, 153, 153, 136, 136, 112, 112, 109, 109, 100, 100, 64,
                                       62, 23, 23, 5, 5, 2, 2))
  plot(evoPhylo,show.node.label = TRUE, cex = 0.6, root.edge = TRUE, tip.color = 'blue')


  # test multiple trees
  two_tree <- data.frame(at = c(0,1,2,2,3,4,0,1,2,3),
                         sus = c(2,3,4,5,6,7,9,10,11,12),
                         inf = c(1,1,3,4,3,5,8,9,9,10))
  class(two_tree) <- c('transmat',class(two_tree))
  expect_message(two_phylo <- as.phylo.transmat(two_tree), regexp = 'found multiple trees')
  expect_true(inherits(two_phylo,'multiPhylo'))
  expect_equal(length(two_phylo),2)
  expect_equal(class(two_phylo[[1]]),'phylo')
  expect_equal(names(two_phylo),c("seed_1", "seed_8"))
  plot(two_phylo[[1]])
  plot(two_phylo[[2]])

  two_phylo_exits <- as.phylo.transmat(two_tree,
                                       vertex.exit.times = c(NA,NA,5,7,NA,NA,5,3,NA,NA,3,NA))
})

test_that("as.network.transmat", {
  skip_on_cran()
  net <- as.network(tm)
  expect_equal(network.size(net), 54)
  expect_equal(list.edge.attributes(net),
               c("actRate", "at", "finalProb", "infDur", "na", "transProb"))
  expect_equal(list.vertex.attributes(net),c("at", "na", "vertex.names"))
})

test_that("plot.transmat", {
  skip_on_cran()
  plot(tm) #phylogram default
  plot(tm, style = "network")
  # plot(tm, style = "gv_tree")
  plot(tm, style = "transmissionTimeline")
})

test_that("'tergmLite == TRUE' correctly errors", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, verbose = FALSE,
                         tergmLite = TRUE)
  mod2 <- netsim(est1, param, init, control)
  expect_error(get_transmat(mod2),"transmat not saved when 'tergmLite == TRUE', check control.net settings")
})

