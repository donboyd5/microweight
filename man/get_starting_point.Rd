% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_poisson.r
\name{get_starting_point}
\alias{get_starting_point}
\title{Construct a starting point from a sample of the data}
\usage{
get_starting_point(p)
}
\arguments{
\item{p}{a problem}
}
\value{
list
}
\description{
\code{get_starting_point} returns a list with a starting point constructed
from a sample and also returns the optimization results from that sample.
}
\examples{
p <- make_problem(h=100, s=3, k=2)
spoint <- get_starting_point(p)
}
