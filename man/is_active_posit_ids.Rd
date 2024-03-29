% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.accessor.R
\name{is_active_posit_ids}
\alias{is_active_posit_ids}
\title{Are These Nodes Active (Positional IDs)}
\usage{
is_active_posit_ids(dat, posit_ids)
}
\arguments{
\item{dat}{Main \code{netsim_dat} object containing a \code{networkDynamic}
object and other initialization information passed from
\code{\link{netsim}}.}

\item{posit_ids}{A vector of node positional identifiers.}
}
\value{
A logical vector with TRUE if the node is still active and FALSE
otherwise.
}
\description{
Are These Nodes Active (Positional IDs)
}
