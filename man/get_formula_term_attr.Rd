% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.utils.R
\name{get_formula_term_attr}
\alias{get_formula_term_attr}
\title{Output ERGM Formula Attributes into a Character Vector}
\usage{
get_formula_term_attr(form, nw)
}
\arguments{
\item{form}{An ERGM model formula.}

\item{nw}{A network object.}
}
\value{
A character vector of vertex attributes.
}
\description{
Given a formation formula for a network model, outputs a
character vector of vertex attributes to be used in
\code{\link{netsim}} simulations.
}
