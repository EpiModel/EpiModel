context("transmat functions")

# test that transmat output has right class
library(EpiModel)
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 0.5)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose.int = 0, use.pids = FALSE)
mod1 <- netsim(est1, param, init, control)
tm <- get_transmat(mod1)

test_that('transmat class',{
  expect_true(is.transmat(tm))
  expect_false(is.transmat(list()))
  expect_equal(class(tm),c('transmat','data.frame'))
  expect_equal(names(tm),c( "at","sus","inf","infDur","transProb", "actRate","finalProb"))
})

# define a known transmat output
tm<-structure(list(at = c(2L, 2L, 4L, 6L, 10L, 11L, 11L, 12L, 12L,  13L, 13L, 13L, 13L, 14L, 14L, 15L, 15L, 16L, 17L, 19L, 19L, 21L,  22L, 22L, 23L, 24L, 25L, 25L, 25L, 26L, 26L, 26L, 26L, 27L, 27L,  31L, 31L, 31L, 32L, 32L, 32L, 32L, 35L, 35L, 36L, 36L, 37L, 37L,  37L, 39L, 40L, 40L, 40L), sus = c(83L, 24L, 75L, 5L, 76L, 9L,  8L, 10L, 38L, 36L, 35L, 21L, 97L, 87L, 66L, 39L, 17L, 78L, 98L,  92L, 82L, 70L, 81L, 60L, 15L, 18L, 51L, 90L, 28L, 14L, 45L, 99L,  56L, 7L, 73L, 27L, 49L, 64L, 46L, 62L, 68L, 42L, 44L, 79L, 69L,  63L, 86L, 25L, 1L, 48L, 32L, 3L, 37L), inf = c(72L, 72L, 72L,  83L, 83L, 83L, 76L, 8L, 8L, 10L, 10L, 10L, 75L, 8L, 36L, 10L,  36L, 66L, 17L, 78L, 36L, 92L, 83L, 70L, 60L, 24L, 39L, 70L, 18L,  70L, 51L, 81L, 51L, 18L, 36L, 36L, 81L, 8L, 49L, 21L, 14L, 8L,  62L, 17L, 38L, 44L, 69L, 73L, 69L, 99L, 35L, 18L, 1L), infDur = c(29,  29, 31, 4, 8, 9, 1, 1, 1, 1, 1, 1, 9, 3, 1, 3, 2, 2, 2, 3, 6,  2, 20, 1, 1, 22, 10, 4, 1, 5, 1, 4, 1, 3, 14, 18, 9, 20, 1, 19,  6, 21, 3, 20, 24, 1, 1, 10, 1, 13, 27, 16, 3), transProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ), actRate = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), finalProb = c(0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 )), .Names = c("at", "sus", "inf", "infDur", "transProb", "actRate",  "finalProb"), row.names = c(NA, 53L), class = c("transmat", "data.frame" ))

# check that it reconstructs tm into the expected phylo object
test_that('as.phylo.transmat',{
  phy<-as.phylo.transmat(tm)
  expect_equal(class(phy),'phylo')
  expect_equal(names(phy),c("edge","Nnode","tip.label","node.label","root.edge","edge.length"))
  # lazy test, this is just veryify that the output matches previous output.
  expect_equal(phy,structure(list(edge = structure(c(26, 26, 26, 27, 27, 27, 30,  30, 31, 31, 31, 30, 32, 31, 32, 34, 32, 27, 36, 36, 29, 36, 35,  28, 35, 29, 32, 32, 28, 30, 30, 34, 33, 33, 29, 27, 29, 3, 1,  30, 2, 31, 33, 32, 23, 19, 4, 36, 35, 34, 5, 6, 28, 7, 8, 9,  16, 10, 22, 11, 12, 21, 13, 15, 14, 17, 18, 20, 25, 24), .Dim = c(35L,  2L)), Nnode = 11L, tip.label = c(5L, 9L, 97L, 87L, 98L, 82L,  15L, 90L, 28L, 45L, 56L, 7L, 27L, 64L, 46L, 68L, 42L, 79L, 63L,  86L, 25L, 48L, 32L, 3L, 37L), node.label = c(72L, 83L, 81L, 18L,  8L, 10L, 36L, 69L, 17L, 51L, 70L), root.edge = 0, edge.length = c(2L,  24L, 13L, 4L, 9L, 9L, 1L, 25L, 1L, 28L, 24L, 3L, 8L, 13L, 2L,  2L, 6L, 20L, 2L, 4L, 1L, 11L, 1L, 17L, 1L, 3L, 24L, 18L, 10L,  20L, 21L, 20L, 1L, 4L, 16L)), .Names = c("edge", "Nnode", "tip.label",  "node.label", "root.edge", "edge.length"), class = "phylo"))
})


test_that('as.network.transmat',{
  net<-as.network.transmat(tm)
  expect_equal(network.size(net),54)
  expect_equal(list.edge.attributes(net),c("actRate","at","finalProb", "infDur","na","transProb"))
  expect_equal(list.vertex.attributes(net),c("at","na","vertex.names"))
})

test_that('plot.transmat',{
  plot(tm) #phylogram default
  plot(tm,style = 'network')
  plot(tm,style='gv_tree')
  plot(tm,style='transmissionTimeline')
})

