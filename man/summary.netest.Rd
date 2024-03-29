% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.netest}
\alias{summary.netest}
\title{Summary for Network Model Fit}
\usage{
\method{summary}{netest}(object, ...)
}
\arguments{
\item{object}{An \code{EpiModel} object of class \code{netest}.}

\item{...}{Additional summary function arguments.}
}
\description{
Prints the summary model fit statistics for an ERGM or STERGM
fit.
}
\details{
This function is simply a wrapper function for \code{summary.ergm}.
Additionally, if the edges dissolution approximation was used to fit the
temporal ERGM, then the dissolution coefficient information will be printed.

If the \code{fit} object is attached to the \code{netest} object, then
\code{summary.netest} will call \code{summary} on \code{fit} using the
\code{...} passed to \code{summary.netest}.  Otherwise,
\code{summary.netest} will print the stored summary of the fit generated
in the original \code{netest} call, using the \code{...} passed to
\code{netest}.
}
\keyword{internal}
