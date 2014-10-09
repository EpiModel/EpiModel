
########################################################################
##
## Testing EpiModel output
##
##  The main testing functions below check for consistency of output:
##  that the compartment sizes and flows are consistent according to
##  discrete-time difference equations
##
########################################################################


#' @title Write Out Test Progress to Console
#'
#' @description Writes the name of a test and \code{...} to console for
#'              showing testing progress.
#'
#' @param test a character string with the name of a test.
#'
#' @keywords internal
#' @export
#'
mcat <- function(test) cat("\n", test, " ... ", sep="")


#' @title Test the Model Output from a Network Model
#'
#' @description Tests whether the model output from a network model is consistent
#'              with key balancing equations for compartment and flow sizes for
#'              each simulation for each time step.
#'
#' @param x an object of class \code{netsim}.
#'
#' @keywords internal
#' @export
#'
test_net <- function(x) {

  type <- x$control$type
  nsims <- x$control$nsims
  vital <- x$param$vital
  modes <- x$param$modes

  for (s in 1:nsims) {
    df <- as.data.frame(x, out = "vals", sim = s)
    for (i in 2:nrow(df)) {

      if (type == "SI") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (modes == 2) {
            test <- with(df, s.num.m2[i] == s.num.m2[i-1] - si.flow.m2[i])
            if (test == FALSE) cat("\nFailed *s.num.m2* at SIM", s, "TIME", i)
            test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i])
            if (test == FALSE) cat("\nFailed *i.num.m2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + b.flow[i] - ds.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (modes == 2) {
            test <- df$num.m2[i] == df$num.m2[i-1] + df$b.flow.m2[i] -
              df$ds.flow.m2[i] - df$di.flow.m2[i]
            if (test == FALSE) cat("\nFailed *num.m2* at SIM", s, "TIME", i)
            test <- df$s.num.m2[i] == df$s.num.m2[i-1] - df$si.flow.m2[i] +
              df$b.flow.m2[i] - df$ds.flow.m2[i]
            if (test == FALSE) cat("\nFailed *s.num.m2* at SIM", s, "TIME", i)
            test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i] - di.flow.m2[i])
            if (test == FALSE) cat("\nFailed *i.num.m2* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIR") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - ir.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i-1] + ir.flow[i])
          if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
          if (modes == 2) {
            test <- with(df, s.num.m2[i] == s.num.m2[i-1] - si.flow.m2[i])
            if (test == FALSE) cat("\nFailed *s.num.m2* at SIM", s, "TIME", i)
            test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i] - ir.flow.m2[i])
            if (test == FALSE) cat("\nFailed *i.num.m2* at SIM", s, "TIME", i)
            test <- with(df, r.num.m2[i] == r.num.m2[i-1] + ir.flow.m2[i])
            if (test == FALSE) cat("\nFailed *r.num.m2* at SIM", s, "TIME", i)
          }
          if (vital == TRUE) {
            test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] -
                           di.flow[i] - dr.flow[i])
            if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + b.flow[i] - ds.flow[i])
            if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - ir.flow[i] - di.flow[i])
            if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
            test <- with(df, r.num[i] == r.num[i-1] + ir.flow[i] - dr.flow[i])
            if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
            if (modes == 2) {
              test <- with(df, num.m2[i] == num.m2[i-1] + b.flow.m2[i] - ds.flow.m2[i] -
                             di.flow.m2[i] - dr.flow.m2[i])
              if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
              test <- with(df, s.num.m2[i] == s.num.m2[i-1] - si.flow.m2[i] +
                             b.flow.m2[i] - ds.flow.m2[i])
              if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
              test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i] -
                             ir.flow.m2[i] - di.flow.m2[i])
              if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
              test <- with(df, r.num.m2[i] == r.num.m2[i-1] + ir.flow.m2[i] - dr.flow.m2[i])
              if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
            }
          }
        }
      }

      if (type == "SIS") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + is.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - is.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (modes == 2) {
            test <- with(df, s.num.m2[i] == s.num.m2[i-1] - si.flow.m2[i] + is.flow.m2[i])
            if (test == FALSE) cat("\nFailed *s.num.m2* at SIM", s, "TIME", i)
            test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i] - is.flow.m2[i])
            if (test == FALSE) cat("\nFailed *i.num.m2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + is.flow[i] +
                         b.flow[i] - ds.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - is.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (modes == 2) {
            test <- with(df, num.m2[i] == num.m2[i-1] + b.flow.m2[i] -
                           ds.flow.m2[i] - di.flow.m2[i])
            if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num.m2[i] == s.num.m2[i-1] - si.flow.m2[i] +
                           is.flow.m2[i] + b.flow.m2[i] - ds.flow.m2[i])
            if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num.m2[i] == i.num.m2[i-1] + si.flow.m2[i] -
                           is.flow.m2[i] - di.flow.m2[i])
            if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          }
        }
      }

    }
  }

}


#' @title Test the Model Output from a Stochastic Individual Contact Model
#'
#' @description Tests whether the model output from an individual contact model
#'              is consistent with key balancing equations for compartment and
#'              flow sizes for each simulation for each time step.
#'
#' @param x an object of class \code{icm}.
#'
#' @keywords internal
#' @export
#'
test_icm <- function(x) {

  type <- x$control$type
  nsims <- x$control$nsims
  vital <- x$param$vital
  groups <- x$param$groups

  for (s in 1:nsims) {
    df <- as.data.frame(x, out = "vals", sim = s)
    for (i in 2:nrow(df)) {

      if (type == "SI") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i-1] - si.flow.g2[i])
            if (test == FALSE) cat("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i])
            if (test == FALSE) cat("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + b.flow[i] - ds.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- df$num.g2[i] == df$num.g2[i-1] + df$b.flow.g2[i] -
              df$ds.flow.g2[i] - df$di.flow.g2[i]
            if (test == FALSE) cat("\nFailed *num.g2* at SIM", s, "TIME", i)
            test <- df$s.num.g2[i] == df$s.num.g2[i-1] - df$si.flow.g2[i] +
              df$b.flow.g2[i] - df$ds.flow.g2[i]
            if (test == FALSE) cat("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) cat("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIR") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - ir.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i-1] + ir.flow[i])
          if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i-1] - si.flow.g2[i])
            if (test == FALSE) cat("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i] - ir.flow.g2[i])
            if (test == FALSE) cat("\nFailed *i.num.g2* at SIM", s, "TIME", i)
            test <- with(df, r.num.g2[i] == r.num.g2[i-1] + ir.flow.g2[i])
            if (test == FALSE) cat("\nFailed *r.num.g2* at SIM", s, "TIME", i)
          }
          if (vital == TRUE) {
            test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] -
                           di.flow[i] - dr.flow[i])
            if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + b.flow[i] - ds.flow[i])
            if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - ir.flow[i] - di.flow[i])
            if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
            test <- with(df, r.num[i] == r.num[i-1] + ir.flow[i] - dr.flow[i])
            if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
            if (groups == 2) {
              test <- with(df, num.g2[i] == num.g2[i-1] + b.flow.g2[i] - ds.flow.g2[i] -
                             di.flow.g2[i] - dr.flow.g2[i])
              if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
              test <- with(df, s.num.g2[i] == s.num.g2[i-1] - si.flow.g2[i] +
                             b.flow.g2[i] - ds.flow.g2[i])
              if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
              test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i] -
                             ir.flow.g2[i] - di.flow.g2[i])
              if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
              test <- with(df, r.num.g2[i] == r.num.g2[i-1] + ir.flow.g2[i] - dr.flow.g2[i])
              if (test == FALSE) cat("\nFailed *r.num* at SIM", s, "TIME", i)
            }
          }
        }
      }

      if (type == "SIS") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + is.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - is.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i-1] - si.flow.g2[i] + is.flow.g2[i])
            if (test == FALSE) cat("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i] - is.flow.g2[i])
            if (test == FALSE) cat("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i-1] + b.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i-1] - si.flow[i] + is.flow[i] +
                         b.flow[i] - ds.flow[i])
          if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i-1] + si.flow[i] - is.flow[i] - di.flow[i])
          if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, num.g2[i] == num.g2[i-1] + b.flow.g2[i] -
                           ds.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) cat("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num.g2[i] == s.num.g2[i-1] - si.flow.g2[i] +
                           is.flow.g2[i] + b.flow.g2[i] - ds.flow.g2[i])
            if (test == FALSE) cat("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i-1] + si.flow.g2[i] -
                           is.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) cat("\nFailed *i.num* at SIM", s, "TIME", i)
          }
        }
      }

    }
  }

}
