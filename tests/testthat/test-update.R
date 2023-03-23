context("update.R functionality")

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
