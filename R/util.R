
#' @title networkLite Constructor Utilities
#'
#' @description Constructor methods for \code{networkLite} objects.
#'
#' @param x either an \code{edgelist} class network representation (including network
#'        attributes in its \code{attributes} list), or a number specifying the network
#'        size.
#' @param attr a named list of vertex attributes for the network represented by \code{x}.
#' @param directed,bipartite,loops,hyper,multiple common network attributes that
#'        may be set via arguments to the \code{networkLite.numeric} method.
#' @param ... additional arguments used by other methods.
#'
#' @details Currently there are two distinct \code{networkLite} constructor methods available.
#'
#'          The \code{edgelist} method takes an \code{edgelist} class object \code{x} with network
#'          attributes attached in its \code{attributes} list, and a named list of vertex
#'          attributes \code{attr}, and returns a \code{networkLite} object, which is a named
#'          list with fields \code{el}, \code{attr}, and \code{gal}; the fields \code{el} and \code{attr} match
#'          the arguments \code{x} and \code{attr} respectively, and the field \code{gal} is the list
#'          of network attributes (copied from \code{attributes(x)}). Missing attributes
#'          \code{directed}, \code{bipartite}, \code{loops}, \code{hyper}, and \code{multiple} are defaulted to
#'          \code{FALSE}; the network size attribute \code{n} must not be missing.  Attributes
#'          \code{class}, \code{dim}, and \code{vnames} (if present) are not copied from \code{x} to the
#'          \code{networkLite}.  (For convenience, a \code{matrix} method, identical to the
#'          \code{edgelist} method, is also defined, to handle cases where the edgelist is,
#'          for whatever reason, not classed as an \code{edgelist}.)
#'
#'          The \code{numeric} method takes a number \code{x} as well as the network attributes
#'          \code{directed}, \code{bipartite}, \code{loops}, \code{hyper}, and \code{multiple} (defaulting to
#'          \code{FALSE}), and returns an empty \code{networkLite} with these network attributes
#'          and number of nodes \code{x}.
#'
#'          Within \code{tergmLite}, the \code{networkLite} data structure is used in the
#'          calls to \code{ergm} and \code{tergm} \code{simulate} functions.
#' 
#' @return  A networkLite object with edge list \code{el}, vertex attributes \code{attr}, and
#'          network attributes \code{gal}.
#'
#' @rdname networkLite
#' @export
#'
#' @examples
#' \dontrun{
#' library("EpiModel")
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # networkLite representation after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#'
#' # Conversion to networkLite class format
#' nwl <- networkLite(dat$el[[1]], dat$attr)
#' nwl
#' }
#'
networkLite <- function(x, ...) {
  UseMethod("networkLite")
}

#' @rdname networkLite
#' @export
networkLite.numeric <- function(x, directed = FALSE, bipartite = FALSE, loops = FALSE, hyper = FALSE, multiple = FALSE, ...) {
  x <- as.numeric(x) # so it's not of class integer

  el <- matrix(0L, nrow = 0L, ncol = 2L)
  attr <- list()
  gal <- list(n = x, directed = directed, bipartite = bipartite, loops = loops, hyper = hyper, multiple = multiple)

  nw <- list(el = el, attr = attr, gal = gal)

  class(nw) <- c("networkLite", "network")
  return(nw)
}

#' @rdname networkLite
#' @export
networkLite.edgelist <- function(x, attr = list(), ...) {
  nw <- list(el = x,
             attr = attr,
             gal = attributes(x)[setdiff(names(attributes(x)), c("class", "dim", "vnames"))])

  # network size attribute is required
  if (is.null(nw$gal[["n"]])) {
    stop("edgelist passed to networkLite constructor must have the `n` attribute.")
  }
  # other common attributes default to FALSE
  if (is.null(nw$gal[["directed"]])) {
    nw$gal[["directed"]] <- FALSE
  }
  if (is.null(nw$gal[["bipartite"]])) {
    nw$gal[["bipartite"]] <- FALSE
  }
  if (is.null(nw$gal[["loops"]])) {
    nw$gal[["loops"]] <- FALSE
  }
  if (is.null(nw$gal[["hyper"]])) {
    nw$gal[["hyper"]] <- FALSE
  }
  if (is.null(nw$gal[["multiple"]])) {
    nw$gal[["multiple"]] <- FALSE
  }

  ## for consistency with network,
  ## we want nw$gal[["n"]] to be of
  ## type numeric, not integer
  nw$gal[["n"]] <- as.numeric(nw$gal[["n"]])

  ## some "double" edgelists have been coming through...
  storage.mode(nw$el) <- "integer"

  class(nw) <- c("networkLite", "network")
  return(nw)
}

#' @rdname networkLite
#' @export
networkLite.matrix <- networkLite.edgelist

#' @name networkLitemethods
#' @title networkLite Methods
#'
#' @description S3 methods for networkLite class, for generics defined in
#'              network package.
#'
#' @param x a \code{networkLite} object.
#' @param attrname the name of an attribute in \code{x}.
#' @param ... any additional arguments.
#'
#' @details Allows use of networkLite objects in \code{ergm_model}.
#'
#' @rdname networkLitemethods
#' @export
#'
get.vertex.attribute.networkLite <- function(x, attrname, ...) {
  x$attr[[attrname]]
}

#' @rdname networkLitemethods
#' @param value attribute value.
#' @param v indices at which to set vertex attribute values.
#' @export
#'
set.vertex.attribute.networkLite <- function(x, attrname, value, v = seq_len(network.size(x)), ...) {
  xn <- substitute(x)

  if(!(attrname %in% list.vertex.attributes(x))) {
    x$attr[[attrname]] <- rep(NA, length = network.size(x))
  }

  x$attr[[attrname]][v] <- value

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
#'
list.vertex.attributes.networkLite <- function(x, ...) {
  names(x$attr)
}

#' @rdname networkLitemethods
#' @export
#'
get.network.attribute.networkLite <- function(x, attrname, ...) {
  x$gal[[attrname]]
}

#' @rdname networkLitemethods
#' @export
#'
set.network.attribute.networkLite <- function(x, attrname, value, ...) {
  xn <- substitute(x)

  x$gal[[attrname]] <- value

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
#'
list.network.attributes.networkLite <- function(x, ...) {
  names(x$gal)
}

#' @rdname networkLitemethods
#' @export
#'
network.edgecount.networkLite <- function(x, ...) {
  NROW(x$el)
}

#' @importFrom tibble as_tibble
#' @rdname networkLitemethods
#' @param output Type of edgelist to output.
#' @export
#'
as.edgelist.networkLite <- function(x, output = c("matrix", "tibble"), ...) {
  output <- match.arg(output)

  if(output == "matrix") {
    m <- x$el
  } else {
    m <- as_tibble(list(.tail = x$el[,1], .head = x$el[,2]))
    class(m) <- c("tibble_edgelist", "edgelist", class(m))
  }

  attr(m, "n") <- as.integer(network.size(x))
  attr(m, "directed") <- as.logical(is.directed(x))
  bip <- if(is.bipartite(x)) x %n% "bipartite" else FALSE
  attr(m, "bipartite") <- if(is.numeric(bip)) as.integer(bip) else bip
  attr(m, "loops") <- as.logical(has.loops(x))
  attr(m, "vnames") <- network.vertex.names(x)

  return(m)
}

#' @rdname networkLitemethods
#' @param object a \code{networkLite} object
#' @param attr specification of a vertex attribute in \code{object} as
#'             described in \code{\link[ergm]{nodal_attributes}}
#' @export
mixingmatrix.networkLite <- function(object, attr, ...) {
  nw <- object

  all_attr <- ergm_get_vattr(attr, nw, multiple = "paste")

  if(is.bipartite(nw)) {
    row_levels <- sort(unique(ergm_get_vattr(attr, nw, bip = "b1", multiple = "paste")))
    col_levels <- sort(unique(ergm_get_vattr(attr, nw, bip = "b2", multiple = "paste")))
  } else {
    row_levels <- sort(unique(all_attr))
    col_levels <- row_levels
  }

  el <- as.edgelist(nw)

  m <- table(from = factor(all_attr[el[,1]], levels = row_levels),
             to = factor(all_attr[el[,2]], levels = col_levels))

  if(!is.bipartite(nw) && !is.directed(nw)) {
    m <- m + t(m) - diag(diag(m))
  }

  m
}

#' @rdname networkLitemethods
#' @param i,j Nodal indices (must be missing for networkLite method)
#' @param value Value to set edges to (must be FALSE for networkLite method)
#' @export
"[<-.networkLite" <- function(x, i, j, value) {
  if(missing(i) && missing(j) && isTRUE(all(value == FALSE))) {
    x$el <- structure(x$el[NULL,,drop=FALSE], class = class(x$el), n = x$gal$n)
    return(x)
  } else {
    stop("networkLite `[<-` operator only supports removing all edges at this time")
  }
}

#' @rdname networkLitemethods
#' @export
print.networkLite <- function(x, ...) {
  cat("networkLite with properties:\n")
  cat("  Network size:", x$gal$n, "\n")
  cat("  Edge count:", NROW(x$el), "\n")
  cat("  Network attributes:", sort(unique(names(x$gal))), "\n")
  cat("  Vertex attributes:", sort(unique(names(x$attr))), "\n")
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
network.naedgecount.networkLite <- function(x, ...) {
  0 # for now
}

#' @rdname networkLitemethods
#' @param tail vector of tails of edges to add to the networkLite
#' @param head vector of heads of edges to add to the networkLite
#' @param names.eval currently unsupported by add.edges.networkLite
#' @param vals.eval currently unsupported by add.edges.networkLite
#' @param check.unique should a check to ensure uniqueness of edges
#'                     in the final edgelist be performed?
#' @export
add.edges.networkLite <- function(x, tail, head, names.eval = NULL, vals.eval = NULL, ..., check.unique = FALSE) {
  if(!is.null(names.eval) || !is.null(vals.eval)) {
    stop("add.edges.networkLite does not currently support ", sQuote("names.eval"), " or ", sQuote("vals.eval"), " arguments.")
  }

  xn <- substitute(x)

  if(length(tail) > 0) {
    new_el <- rbind(x$el, cbind(tail, head))
    new_el <- new_el[order(new_el[,1], new_el[,2]),]
    if(check.unique) {
      ## this could be made faster by exploiting
      ## the fact that new_el is sorted
      new_el <- unique(new_el)
    }
    class(new_el) <- c("edgelist", class(new_el))
    attr(new_el, "n") <- x$gal$n
    x$el <- new_el
  }

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
as.networkLite <- function(x, ...) {
  UseMethod("as.networkLite")
}

#' @rdname networkLitemethods
#' @export
as.networkLite.network <- function(x, ...) {
  edgelist <- as.edgelist(x)
  vertex_attributes <- list()
  for(name in list.vertex.attributes(x)) {
    vertex_attributes[[name]] <- x %v% name
  }
  rv <- networkLite(edgelist, vertex_attributes)
  for(name in list.network.attributes(x)) {
    rv$gal[[name]] <- x %n% name
  }
  rv
}

#' @rdname networkLitemethods
#' @export
as.networkLite.networkLite <- function(x, ...) {
  x
}
