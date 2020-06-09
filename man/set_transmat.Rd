% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.mod.infection.R
\name{set_transmat}
\alias{set_transmat}
\title{Save Transmission Matrix}
\usage{
set_transmat(dat, del, at)
}
\arguments{
\item{dat}{Master list object containing a \code{networkDynamic} object or
edgelist (if tergmLite is used) and other initialization information
passed from \code{\link{netsim}}.}

\item{del}{Discordant edgelist created within \code{infection.net} and
\code{infection.2g.net}.}

\item{at}{Current time step.}
}
\description{
This function appends the transmission matrix created during
             \code{infection.net} and \code{infection.2g.net}.
}
\details{
This internal function works within the parent \code{\link{infection.net}} functions
to save the transmission matrix created at time step \code{at} to the master list object
\code{dat}.
}