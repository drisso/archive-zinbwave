% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{getTheta}
\alias{getTheta}
\title{Returns the vector of inverse dispersion parameters}
\usage{
getTheta(object)
}
\arguments{
\item{object}{an object that describes a matrix of zero-inflated
distributions.}
}
\value{
the vector of inverse dispersion parameters theta
}
\description{
Given an object that describes a matrix of zero-inflated negative binomial
distributions, returns the vector of inverse dispersion parameters
\code{theta}.
}
\examples{
a <- zinbModel(n=5, J=10)
getTheta(a)
}
