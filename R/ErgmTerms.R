
#' @title Definition for absdiffnodemix ERGM Term
#'
#' @description This function defines and initialize the absdiffnodemix ERGM term
#'              that allows for targeting age homophily by race.
#'
#' @param nw An object of class \code{network}.
#' @param arglist A list of arguments as specified in the \code{ergm.userterms}
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        \code{ergm.userterms} package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in partnership
#' formation that is heterogenous by race. The absdiff component allows targets
#' the distribution of age mixing on that continuous variable, and the nodemix
#' component differentiates this for black-black, black-white, and white-white
#' couples.
#'
#' @author Steven M. Goodreau
#'
#' @export
InitErgmTerm.absdiffnodemix <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attrname", "byattrname"),
                      vartypes = c("character", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  nodecov <- get.node.attr(nw, a$attrname)
  nodecovby <- get.node.attr(nw, a$byattrname)
  u <- sort(unique(nodecovby))
  if (any(is.na(nodecovby))) {
    u <- c(u, NA)
  }

  nodecovby <- match(nodecovby, u, nomatch = length(u) + 1)
  ui <- seq(along = u)

  uui <- matrix(1:length(ui) ^ 2, length(ui), length(ui))
  urm <- t(sapply(ui, rep, length(ui)))
  ucm <- sapply(ui, rep, length(ui))
  uun <- outer(u, u, paste, sep = ".")
  uui <- uui[upper.tri(uui, diag = TRUE)]
  urm <- urm[upper.tri(urm, diag = TRUE)]
  ucm <- ucm[upper.tri(ucm, diag = TRUE)]
  uun <- uun[upper.tri(uun, diag = TRUE)]

  inputs = c(length(nodecov), length(urm), nodecov, nodecovby, urm, ucm)

  list(name = "absdiffnodemix",
       coef.names = paste("absdiffnodemix", a$attrname, a$byattrname, uun, sep = "."),
       pkgname = "EpiModel",
       inputs = inputs,
       dependence = FALSE)
}


#' @title Definition for absdiffby ERGM Term
#'
#' @description This function defines and initialize the absdiffby ERGM term
#'              that allows for targeting age homophily by sex.
#'
#' @param nw An object of class \code{network}.
#' @param arglist A list of arguments as specified in the \code{ergm.userterms}
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        \code{ergm.userterms} package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in partnership
#' formation that is asymetric by sex. The absdiff component targets age homophily
#' while the by component allows that to be structed by a binary attribute such
#' as "male", in order to enforce an offset in the average difference.
#'
#' @author Samuel M. Jenness
#'
#' @export
InitErgmTerm.absdiffby <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attrname", "by", "assym"),
                      vartypes = c("character", "character", "numeric"),
                      required = c(TRUE, TRUE, TRUE),
                      defaultvalues = list(NULL, NULL, NULL))

  nodecov <- get.node.attr(nw, a$attrname)
  nodeby <- get.node.attr(nw, a$by)
  coef.names <- paste("absdiffby", a$attrname, a$by, sep = ".")

  list(name = "absdiffby",
       coef.names = coef.names,
       pkgname = "EpiModel",
       inputs = c(a$assym, nodecov, nodeby),
       dependence = FALSE,
       emptynwstats = 0
  )
}
