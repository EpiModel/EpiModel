
#' @title networkLite Constructor Utilities
#'
#' @description Constructor methods for \code{networkLite} objects.
#'
#' @param x Either an \code{edgelist} class network representation (including
#'        network attributes in its \code{attributes} list), or a number
#'        specifying the network size.
#' @param attr A named list of vertex attributes for the network represented by
#'        \code{x}.
#' @param directed,bipartite Common network attributes that may be set via 
#'        arguments to the \code{networkLite.numeric} method.
#' @param ... Additional arguments used by other methods.
#'
#' @details Currently there are several distinct \code{networkLite} constructor
#' methods available.
#'
#' The \code{edgelist} method takes an \code{edgelist} class object \code{x}
#' with network attributes attached in its \code{attributes} list, and a named
#' list of vertex attributes \code{attr}, and returns a \code{networkLite}
#' object, which is a named list with fields \code{el}, \code{attr}, and
#' \code{gal}; the fields \code{el} and \code{attr} match the arguments \code{x}
#' and \code{attr} (the latter coerced to \code{tibble}) respectively, and the 
#' field \code{gal} is the list of network attributes (copied from 
#' \code{attributes(x)}). Missing network attributes \code{directed} and 
#' \code{bipartite} are defaulted to \code{FALSE}; the network size attribute 
#' \code{n} must not be missing. Attributes \code{class}, \code{dim}, 
#' \code{dimnames}, \code{vnames}, and \code{mnext} (if present) are not copied 
#' from \code{x} to the \code{networkLite}.  (For convenience, a \code{matrix} 
#' method, identical to the \code{edgelist} method, is also defined, to handle 
#' cases where the edgelist is, for whatever reason, not classed as an 
#' \code{edgelist}.)
#'
#' The \code{numeric} method takes a number \code{x} as well as the network
#' attributes \code{directed} and \code{bipartite} (defaulting to \code{FALSE}), 
#' and returns an empty \code{networkLite} with these network attributes and 
#' number of nodes \code{x}.
#'
#' The constructor \code{networkLite_initialize} is also available for creating
#' an empty \code{networkLite}, and its \code{x} argument should be a number
#' indicating the size of the \code{networkLite} to create.
#'
#' Within \code{tergmLite}, the \code{networkLite} data structure is used in the
#' calls to \code{ergm} and \code{tergm} \code{simulate} functions.
#'
#' @return
#' A networkLite object with edge list \code{el}, vertex attributes \code{attr},
#' and network attributes \code{gal}.
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
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        tergmLite = TRUE)
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
networkLite.edgelist <- function(x, 
                                 attr = list(vertex.names = seq_len(attributes(x)[["n"]]), 
                                             na = logical(attributes(x)[["n"]])), 
                                 ...) {
  nw <- list(el = x,
             attr = as_tibble(attr),
             gal = attributes(x)[setdiff(names(attributes(x)),
                                         c("class", "dim", "dimnames", "vnames", "mnext"))])

  if(!is_tibble(x)) {
    nw$el <- as_tibble(list(.tail = as.integer(x[,1]), .head = as.integer(x[,2])))
  }

  nw$el[["na"]] <- NVL(nw$el[["na"]], logical(NROW(nw$el)))
  nw$el[["na"]][is.na(nw$el[["na"]])] <- FALSE

  # network size attribute is required
  if (is.null(nw$gal[["n"]])) {
    stop("edgelist passed to networkLite must have the `n` attribute.")
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

  if (!isFALSE(nw$gal[["loops"]]) || !isFALSE(nw$gal[["hyper"]]) || !isFALSE(nw$gal[["multiple"]])) {
    stop("networkLite requires network attributes `loops`, `hyper`, and `multiple` be `FALSE`.")
  }
  
  ## for consistency with network,
  ## we want nw$gal[["n"]] to be of
  ## type numeric, not integer
  nw$gal[["n"]] <- as.numeric(nw$gal[["n"]])

  class(nw) <- c("networkLite", "network")  
  return(nw)
}

#' @rdname networkLite
#' @export
networkLite.matrix <- networkLite.edgelist

#' @rdname networkLite
#' @export
networkLite.numeric <- function(x,
                                directed = FALSE,
                                bipartite = FALSE,
                                ...) {
  x <- as.numeric(x) # so it's not of class integer

  el <- as_tibble(list(.tail = integer(0), .head = integer(0), na = logical(0)))
  attr <- list(vertex.names = seq_len(x), na = logical(x))
  gal <- list(n = x, directed = directed, bipartite = bipartite,
              loops = FALSE, hyper = FALSE, multiple = FALSE)

  nw <- list(el = el, attr = as_tibble(attr), gal = gal)

  class(nw) <- c("networkLite", "network")
  return(nw)
}

#' @rdname networkLite
#' @export
networkLite_initialize <- networkLite.numeric

#' @name networkLitemethods
#' @title networkLite Methods
#'
#' @description S3 methods for networkLite class, for generics defined in
#'              network package.
#'
#' @param x A \code{networkLite} object.
#' @param attrname The name of an attribute in \code{x}.
#' @param value The attribute value to set in vertex, edge, and network 
#'              attribute setters; the value to set edges to (must be FALSE)
#'              for the \code{networkLite} replacement method.
#' @param ... Any additional arguments.
#'
#' @details Allows use of networkLite objects in \code{ergm_model}.
#'
#' @return An edgelist for \code{as.edgelist.networkLite}; an updated
#'         \code{networkLite} object for the replacement method. The other
#'         methods return no objects.
#'
#' @rdname networkLitemethods
#' @export
#'
get.vertex.attribute.networkLite <- function(x, attrname, ...) {
  x$attr[[attrname]]
}

#' @rdname networkLitemethods
#' @param v Indices at which to set vertex attribute values.
#' @export
#'
set.vertex.attribute.networkLite <- function(x,
                                             attrname,
                                             value,
                                             v = seq_len(network.size(x)),
                                             ...) {
  xn <- substitute(x)

  if (!(attrname %in% list.vertex.attributes(x))) {
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
get.edge.attribute.networkLite <- function(x, attrname, ...) {
  x$el[[attrname]]
}

#' @rdname networkLitemethods
#' @export
#'
get.edge.value.networkLite <- get.edge.attribute.networkLite

#' @rdname networkLitemethods
#' @param e edge indices to assign value
#' @export
#'
set.edge.attribute.networkLite <- function(x, attrname, value, e = seq_len(network.edgecount(x, na.omit = FALSE)), ...) {
  xn <- substitute(x)

  if (!(attrname %in% list.edge.attributes(x))) {
    x$el[[attrname]] <- rep(NA, length = network.edgecount(x, na.omit = FALSE))
  }

  x$el[[attrname]][e] <- value

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
#'
set.edge.value.networkLite <- function(x, attrname, value, e = seq_len(network.edgecount(x, na.omit = FALSE)), ...) {
  xn <- substitute(x)

  if (!(attrname %in% list.edge.attributes(x))) {
    x$el[[attrname]] <- rep(NA, length = network.edgecount(x, na.omit = FALSE))
  }

  x$el[[attrname]][e] <- value[as.matrix(x$el[e,c(".tail", ".head")])]

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)  
}

#' @rdname networkLitemethods
#' @export
#'
list.edge.attributes.networkLite <- function(x, ...) {
  colnames(x$el)[-c(1,2)]
}


#' @rdname networkLitemethods
#' @param na.omit logical; omit missing edges from edge count?
#' @export
#'
network.edgecount.networkLite <- function(x, na.omit = TRUE, ...) {
  if (na.omit == TRUE) {
    NROW(x$el) - network.naedgecount(x)
  } else {
    NROW(x$el)
  }
}

#' @rdname networkLitemethods
#' @param output Type of edgelist to output.
#' @param na.rm should missing edges be dropped from edgelist?
#' @export
#'
as.edgelist.networkLite <- function(x, attrname = NULL, output = c("matrix", "tibble"), na.rm = TRUE, ...) {
  output <- match.arg(output)

  if (output == "matrix") {
    m <- matrix(c(x$el$.tail, x$el$.head), ncol = 2)
    if(!is.null(attrname)) {
      m <- cbind(m, x$el[[attrname]])
    }
  } else {
    m <- x$el[c(".tail", ".head", attrname)]
  }

  if (na.rm && NROW(m) > 0) {
    na <- NVL(x %e% "na", FALSE)
    m <- m[!na,,drop=FALSE]
  }

  attr(m, "dimnames") <- NULL
  
  attr(m, "n") <- as.integer(network.size(x))
  attr(m, "vnames") <- network.vertex.names(x)
  bip <- if (is.bipartite(x)) x %n% "bipartite" else FALSE
  attr(m, "bipartite") <- if (is.numeric(bip)) as.integer(bip) else bip
  attr(m, "directed") <- as.logical(is.directed(x))
  attr(m, "loops") <- as.logical(has.loops(x))
  class(m) <- c(if(output == "matrix") "matrix_edgelist" else "tibble_edgelist", "edgelist", class(m))
  return(m)
}

#' @rdname networkLitemethods
#' @param object A \code{networkLite} object.
#' @param attr Specification of a vertex attribute in \code{object} as
#'             described in \code{\link[ergm]{nodal_attributes}}.
#' @export
mixingmatrix.networkLite <- function(object, attr, ...) {  
  nw <- object

  all_attr <- ergm_get_vattr(attr, nw, multiple = "paste")

  if (is.bipartite(nw)) {
    row_levels <- sort(unique(ergm_get_vattr(attr, nw, bip = "b1",
                                             multiple = "paste")))
    col_levels <- sort(unique(ergm_get_vattr(attr, nw, bip = "b2",
                                             multiple = "paste")))
  } else {
    row_levels <- sort(unique(all_attr))
    col_levels <- row_levels
  }

  el <- as.edgelist(nw)

  m <- table(from = factor(all_attr[el[, 1]], levels = row_levels),
             to = factor(all_attr[el[, 2]], levels = col_levels))

  if (!is.bipartite(nw) && !is.directed(nw)) {
    m <- m + t(m) - diag(diag(m))
  }

  m
}

#' @rdname networkLitemethods
#' @param i,j Nodal indices (must be missing for networkLite method).
#' @param add.edges should edges being assigned to be added if not already 
#'        present?
#' @export
"[<-.networkLite" <- function(x, i, j, names.eval = NULL, add.edges = FALSE, value) {
  if (!missing(i) || !missing(j)) {
    stop("`[<-.networkLite` does not support `i` and `j` arguments at this time")
  }
  
  if (any(is.na(value))) {
    stop("`[<-.networkLite` does not support NA `value` arguments at this time")
  }
  
  if (is.null(names.eval) && isTRUE(all(value == FALSE))) {
    x$el <- as_tibble(list(.tail = integer(0), .head = integer(0), na = logical(0)))
    return(x)
  }
  
  b1 <- if (is.bipartite(x)) x %n% "bipartite" else network.size(x)
  b2 <- if (is.bipartite(x)) network.size(x) - x %n% "bipartite" else network.size(x)
  
  if (!is.matrix(value)) {
    value <- matrix(value, nrow = b1, ncol = b2)
  } else {
    if (nrow(value) < b1 || ncol(value) < b2) {
      stop("too small a matrix `value` passed to `[<-.networkLite`")
    }
    value <- value[seq_len(b1),seq_len(b2),drop=FALSE]
  }
  
  if (is.null(names.eval)) {
    # add edges whether or not add.edges is TRUE, for consistency with `network` behavior
    w <- which(value != 0, arr.ind = TRUE)
    if (is.bipartite(x)) {
      w[,2] <- w[,2] + b1
    }
    if (!is.directed(x)) {
      w <- w[w[,1] < w[,2],,drop=FALSE]
    } else {
      w <- w[w[,1] != w[,2],,drop=FALSE]
    }
    w <- w[order(w[,1], w[,2]),,drop=FALSE]
    x$el <- as_tibble(list(.tail = w[,1], .head = w[,2], na = logical(NROW(w))))
  } else {
    if (!add.edges) {
      el <- as.edgelist(x, na.rm = FALSE)
      if (is.bipartite(x)) {
        el[,2] <- el[,2] - b1
      }
      if (names.eval == "na") {
        value[is.na(value)] <- FALSE
      }
      set.edge.attribute(x, names.eval, value[el])
    } else {
      w <- which(value != 0, arr.ind = TRUE)
      vals <- value[w]
      if (is.bipartite(x)) {
        w[,2] <- w[,2] + b1
      }
      if (!is.directed(x)) {
        vals <- vals[w[,1] < w[,2]]
        w <- w[w[,1] < w[,2],,drop=FALSE]
      } else {
        vals <- vals[w[,1] != w[,2]]
        w <- w[w[,1] != w[,2],,drop=FALSE]
      }
      vals <- vals[order(w[,1], w[,2])]
      w <- w[order(w[,1], w[,2]),,drop=FALSE]
      if (names.eval == "na") {
        vals[is.na(vals)] <- FALSE
        tbl_list <- list(w[,1], w[,2], vals)
        names(tbl_list) <- c(".tail", ".head", names.eval)
      } else {
        tbl_list <- list(w[,1], w[,2], vals, logical(NROW(w)))
        names(tbl_list) <- c(".tail", ".head", names.eval, "na")      
      }
      x$el <- as_tibble(tbl_list)
    }
  }
  return(x)
}

#' @rdname networkLitemethods
#' @export
print.networkLite <- function(x, ...) {
  cat("networkLite with properties:\n")
  cat("  Network size:", network.size(x), "\n")
  cat("  Edge count:", network.edgecount(x, na.omit = FALSE), "\n")  
  cat("    Non-missing edge count:", network.edgecount(x, na.omit = TRUE), "\n")
  cat("    Missing edge count:", network.naedgecount(x), "\n")
  cat("  Network attributes:", list.network.attributes(x), "\n")
  cat("  Vertex attributes:", list.vertex.attributes(x), "\n")
  cat("  Edge attributes:", list.edge.attributes(x), "\n")  
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
network.naedgecount.networkLite <- function(x, ...) {
  sum(x$el[["na"]])
}

#' @rdname networkLitemethods
#' @param tail Vector of tails of edges to add to the networkLite.
#' @param head Vector of heads of edges to add to the networkLite.
#' @param names.eval name(s) of edge attributes
#' @param vals.eval value(s) of edge attributes
#' @export
add.edges.networkLite <- function(x, tail, head, names.eval = NULL,
                                  vals.eval = NULL, ...) {
  tail <- NVL(unlist(tail), integer(0))
  head <- NVL(unlist(head), integer(0))
  if(length(names.eval) == 0) {
    update_tibble <- as_tibble(list(.tail = tail, .head = head, na = logical(length(tail))))
  } else if(is.character(names.eval) && length(names.eval) == 1) {
    tbl_list <- list(tail, head, unlist(vals.eval))
    names(tbl_list) <- c(".tail", ".head", names.eval)
    update_tibble <- as_tibble(tbl_list)
  } else {
    if(!is.list(names.eval)) names.eval <- as.list(rep(names.eval, length = length(tail)))
    if(!is.list(vals.eval)) vals.eval <- as.list(rep(vals.eval, length = length(names.eval)))
    
    for(i in seq_along(vals.eval)) {
      vals.eval[[i]] <- as.list(vals.eval[[i]])
      names(vals.eval[[i]]) <- unlist(names.eval[[i]])
    }
    
    update_tibble <- dplyr::bind_cols(as_tibble(list(.tail = tail, .head = head)),
                                      dplyr::bind_rows(lapply(vals.eval, function(x) if(length(x) > 0) as_tibble(x) else tibble(NULL, .rows = 1))))
  }

  update_tibble[["na"]] <- NVL(update_tibble[["na"]], logical(NROW(update_tibble)))
  update_tibble[["na"]][is.na(update_tibble[["na"]])] <- FALSE

  xn <- substitute(x)

  x$el <- dplyr::bind_rows(x$el, update_tibble)
  x$el <- x$el[order(x$el$.tail, x$el$.head),]
  x$el <- x$el[!duplicated(x$el[,c(".tail", ".head")]),]

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
  el <- as.edgelist(x, na.rm = FALSE)

  rv <- networkLite(el)
  
  for (name in list.vertex.attributes(x)) {
    rv %v% name <- x %v% name
  }
  
  for (name in setdiff(list.network.attributes(x), c("mnext"))) {
    rv %n% name <- x %n% name
  }

  eids <- unlist(get.dyads.eids(x, el[,1], el[,2], na.omit = FALSE))
  for (name in list.edge.attributes(x)) {
    rv %e% name <- unlist(get.edge.attribute(x, name, null.na = TRUE, deleted.edges.omit = FALSE, unlist = FALSE)[eids])
  }
  
  for (name in setdiff(names(attributes(x)), c("class", "names"))) {
    attr(rv, name) <- attr(x, name)
  }
  
  rv
}

#' @rdname networkLitemethods
#' @export
as.networkLite.networkLite <- function(x, ...) {
  x
}

# internal convenience function for converting from networkLite to network
to_network_networkLite <- function(x, ...) {
  nw <- network.initialize(network.size(x), 
                           directed = x %n% "directed",
                           bipartite = x %n% "bipartite")
  
  el <- as.edgelist(x, na.rm = FALSE)
  
  nw <- add.edges(nw, el[,1], el[,2])
  
  for(name in list.vertex.attributes(x)) {
    nw %v% name <- x %v% name
  }
  
  for(name in list.network.attributes(x)) {
    nw %n% name <- x %n% name
  }

  eids <- unlist(get.dyads.eids(nw, el[,1], el[,2], na.omit = FALSE))
  for(name in list.edge.attributes(x)) {
    set.edge.attribute(nw, name, x %e% name, eids)
  }
  
  for (name in setdiff(names(attributes(x)), c("class", "names"))) {
    attr(nw, name) <- attr(x, name)
  }
  
  nw
}

#' @rdname networkLitemethods
#' @export
as.networkDynamic.networkLite <- function(object, ...) {  
  as.networkDynamic(to_network_networkLite(object))
}

#' @rdname networkLitemethods
#' @export
valid.eids.networkLite <- function(x, ...) {
  seq_len(NROW(x$el))
}

#' @rdname networkLitemethods
#' @param attrnames vector specifying edge attributes to include in the tibble;
#'        may be logical, integer, or character vector, the former two being 
#'        used to select attribute names from \code{list.edge.attributes(x)}, 
#'        and the latter being used as the attribute names themselves
#' @export
as_tibble.networkLite <- function(x, attrnames = NULL, na.rm = TRUE, ...) {
  if(is.logical(attrnames) || is.numeric(attrnames)) attrnames <- na.omit(list.edge.attributes(x)[attrnames])
  out <- x$el[,c(".tail", ".head", attrnames)]
  if(na.rm && NROW(out) > 0) {
    na <- NVL(x %e% "na", FALSE)
    out <- out[!na,]
  }
  attr(out, "n") <- network.size(x)
  attr(out, "vnames") <- network.vertex.names(x)
  if(is.bipartite(x)) attr(out, "bipartite") <- x %n% "bipartite"
  out
}

#' @rdname networkLitemethods
#' @param matrix.type type of matrix to return from 
#'        \code{as.matrix.networkLite}
#' @export
as.matrix.networkLite <- function(x, matrix.type = c("adjacency", "incidence", "edgelist"), attrname = NULL, ...) {
  matrix.type <- match.arg(matrix.type)
  switch(matrix.type,
         adjacency = as.matrix.networkLite.adjacency(x, attrname, ...),
         incidence = as.matrix.networkLite.incidence(x, attrname, ...),
         edgelist = as.matrix.networkLite.edgelist(x, attrname, ...))
}

as.matrix.networkLite.adjacency <- function(x, attrname = NULL, ...) {
  el <- as.edgelist(x, na.rm = FALSE)
  
  if(!is.null(attrname)) {
    vals <- x %e% attrname
  } else {
    vals <- rep(1, network.edgecount(x, na.omit = FALSE))
  }
  vals[NVL(x %e% "na", FALSE)] <- NA
  
  n <- network.size(x)

  m <- matrix(0, nrow = n, ncol = n)
  m[el] <- vals
  if(!is.directed(x)) {
    m[el[,c(2,1)]] <- vals
  }
  dimnames(m) <- rep(list(network.vertex.names(x)), 2)
  
  if(is.bipartite(x)) {
    bip <- x %n% "bipartite"
    m[seq_len(bip), -seq_len(bip)]
  } else {
    m
  }
}

as.matrix.networkLite.incidence <- function(x, attrname = NULL, ...) {
  el <- as.edgelist(x, na.rm = FALSE)
  
  vals <- NVL2(attrname, x %e% attrname, rep(1, network.edgecount(x, na.omit = FALSE)))
  vals[NVL(x %e% "na", FALSE)] <- NA
  
  m <- matrix(0, nrow = network.size(x), ncol = network.edgecount(x, na.omit = FALSE))

  m[cbind(el[,1], seq_len(NROW(el)))] <- if(is.directed(x)) -vals else vals
  m[cbind(el[,2], seq_len(NROW(el)))] <- vals

  m
}

as.matrix.networkLite.edgelist <- function(x, attrname = NULL, na.rm = TRUE, ...) {
  m <- matrix(c(x$el$.tail, x$el$.head), ncol = 2)
  if (!is.null(attrname)) {
    m <- cbind(m, x$el[[attrname]])
  }
  if (na.rm == TRUE) {
    m <- m[!NVL(x %e% "na", FALSE),,drop=FALSE]
  }
  attr(m, "n") <- network.size(x)
  attr(m, "vnames") <- network.vertex.names(x)
  if (is.bipartite(x)) attr(m, "bipartite") <- x %n% "bipartite"
  m
}

#' @rdname networkLitemethods
#' @export
is.na.networkLite <- function(x) {
  y <- networkLite(network.size(x),
                   directed = x %n% "directed",
                   bipartite = x %n% "bipartite")
  el <- as.edgelist(x, na.rm = FALSE)
  elna <- el[NVL(x %e% "na", FALSE),,drop=FALSE]
  add.edges(y, elna[,1], elna[,2])
  y
}

#' @rdname networkLitemethods
#' @export
delete.vertex.attribute.networkLite <- function(x, attrname, ...) {
  xn <- substitute(x)

  x$attr[[attrname]] <- NULL

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
delete.edge.attribute.networkLite <- function(x, attrname, ...) {
  xn <- substitute(x)

  x$el[[attrname]] <- NULL

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @export
delete.network.attribute.networkLite <- function(x, attrname, ...) {
  xn <- substitute(x)

  x$gal[[attrname]] <- NULL

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @param eid edge indices to delete (between 1 and 
#'        \code{network.edgecount(x, na.omit = FALSE)})
#' @export
delete.edges.networkLite <- function(x, eid, ...) {
  xn <- substitute(x)

  eid <- as.integer(eid)
  eid <- eid[eid >= 1 & eid <= NROW(x$el)]
  if(length(eid) > 0) {
    x$el <- x$el[-eid,]
  }

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @param vid vertex indices to delete (between 1 and \code{network.size(x)})
#' @export
delete.vertices.networkLite <- function(x, vid, ...) {
  xn <- substitute(x)

  vid <- as.integer(vid)
  vid <- vid[vid >= 1 & vid <= network.size(x)]
  if(length(vid) > 0) {
    # drop edges with deleted nodes
    x$el <- x$el[!(x$el$.tail %in% vid | x$el$.head %in% vid),]

    # drop vertex attributes for deleted nodes
    x$attr <- x$attr[-vid,]

    # remap nodal indices for remaining edges
    a <- seq_len(network.size(x))
    b <- integer(network.size(x))
    b[vid] <- 1L
    b <- cumsum(b)
    a <- a - b
    x$el$.tail <- a[x$el$.tail]
    x$el$.head <- a[x$el$.head]

    # update network attributes
    x %n% "n" <- x %n% "n" - length(vid)
    if(is.bipartite(x)) {
      x %n% "bipartite" <- x %n% "bipartite" - sum(vid <= x %n% "bipartite")
    }    
  }

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @param nv number of vertices to add to the \code{networkLite}
#' @param vattr list (of length \code{nv}) of named lists of vertex attributes 
#'        for added vertices, or \code{NULL} to indicate vertex attributes are
#'        not being passed
#' @param last.mode logical; if \code{x} is bipartite, should the new vertices
#'        be added to the second mode?
#' @export
add.vertices.networkLite <- function(x, nv, vattr = NULL, last.mode = TRUE, ...) {
  xn <- substitute(x)

  nv <- as.integer(nv)
  if(nv > 0) {
    oldsize <- network.size(x)
    x %n% "n" <- oldsize + nv

    if(is.bipartite(x) && !last.mode) {
      offset <- x %n% "bipartite"
      x %n% "bipartite" <- x %n% "bipartite" + nv
      x$el$.head <- x$el$.head + nv
    } else {
      offset <- oldsize
    }

    if(!is.null(vattr)) {
      update_tibble <- dplyr::bind_rows(lapply(vattr, function(x) if(length(x) > 0) as_tibble(x) else tibble(NULL, .rows=1)))
    } else {
      update_tibble <- as_tibble(list(na = logical(nv)))
    }
    update_tibble[["na"]] <- NVL(update_tibble[["na"]], logical(NROW(update_tibble)))
    update_tibble[["na"]][is.na(update_tibble[["na"]])] <- FALSE
    
    x$attr <- dplyr::bind_rows(x$attr[seq_len(offset),],
                               update_tibble,
                               x$attr[offset + seq_len(oldsize - offset),])
  }

  on.exit(eval.parent(call("<-", xn, x)))
  invisible(x)
}

#' @rdname networkLitemethods
#' @param e1,e2 networkLite objects
#' @export
`+.networkLite` <- function(e1, e2) {
  if (!identical(e1 %n% "n", e2 %n% "n") ||
      !identical(e1 %n% "directed", e2 %n% "directed") ||
      !identical(e1 %n% "bipartite", e2 %n% "bipartite")) {
    stop("cannot add networkLites of differing network size, directedness, or bipartiteness")
  }

  if (any(NVL(e1 %e% "na", FALSE)) || any(NVL(e2 %e% "na", FALSE))) {
    stop("adding networkLites with missing edges is not currently supported")
  }
  
  out <- e1
  if(network.edgecount(e2, na.omit = FALSE) > 0) {
    edgelist <- dplyr::bind_rows(e1$el, e2$el)
    edgelist <- edgelist[!duplicated(edgelist[,c(".tail", ".head")]),]
    out$el <- edgelist[order(edgelist[,1], edgelist[,2]),]
  }
  out
}

#' @rdname networkLitemethods
#' @export
`-.networkLite` <- function(e1, e2) {
  if (!identical(e1 %n% "n", e2 %n% "n") ||
      !identical(e1 %n% "directed", e2 %n% "directed") ||
      !identical(e1 %n% "bipartite", e2 %n% "bipartite")) {
    stop("cannot subtract networkLites of differing network size, directedness, or bipartiteness")
  }

  if (any(NVL(e1 %e% "na", FALSE)) || any(NVL(e2 %e% "na", FALSE))) {
    stop("subtracting networkLites with missing edges is not currently supported")
  }

  out <- e1
  if(network.edgecount(e2, na.omit = FALSE) > 0) {
    edgelist <- dplyr::bind_rows(e2$el, e1$el)
    nd <- !duplicated(edgelist[,c(".tail", ".head")])
    out$el <- out$el[nd[-seq_len(network.edgecount(e2, na.omit = FALSE))],]
    out$el <- out$el[order(out$el[,1], out$el[,2]),]
  }
  out
}
