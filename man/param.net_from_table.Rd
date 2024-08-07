% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_df.R
\name{param.net_from_table}
\alias{param.net_from_table}
\title{Parameters List for Stochastic Network Models from a Formatted
Data Frame}
\usage{
param.net_from_table(long.param.df)
}
\arguments{
\item{long.param.df}{A \code{data.frame} of parameters. See details for the
expected format.}
}
\value{
A list object of class \code{param.net}, which can be passed to
\code{\link{netsim}}.
}
\description{
Sets the epidemic parameters for stochastic network models with
\code{\link{netsim}} using a specially formatted data frame of
parameters.
}
\section{long.param.df}{

It is possible to set input parameters using a specifically formatted
\code{data.frame} object. The first 3 columns of this \code{data.frame} must
be:
\itemize{
\item \code{param}: The name of the parameter. If this is a non-scalar
parameter (a vector of length > 1), end the parameter name with the
position on the vector (e.g., \code{"p_1"}, \code{"p_2"}, ...).
\item \code{value}: the value for the parameter (or the value of the
parameter in the Nth position if non-scalar).
\item \code{type}: a character string containing either \code{"numeric"},
\code{"logical"}, or \code{"character"} to define the parameter object
class.
}

In addition to these 3 columns, the \code{data.frame} can contain any number
of other columns, such as \code{details} or \code{source} columns to document
parameter meta-data. However, these extra columns will not be used by
EpiModel.
}

