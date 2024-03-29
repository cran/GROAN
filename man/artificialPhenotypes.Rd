% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/artificial_data.R
\name{artificialPhenotypes}
\alias{artificialPhenotypes}
\title{Create a set of artificial phenotypes}
\usage{
artificialPhenotypes(
  genotypes,
  mu = 0,
  h2 = 0.8,
  markersEffectDistr = "rnorm",
  ...
)
}
\arguments{
\item{genotypes}{a samples x markers matrix of 0/1/2...ploidy}

\item{mu}{intercept (added to everything)}

\item{h2}{heritability (must be greater than zero, an less than or equal to one)}

\item{markersEffectDistr}{name of the random generation function of
the selected distribution for markers effects. It must
accept n as argument indicating the number of elements
to be returned.}

\item{...}{further parameters are passed to markersEffectDistr}
}
\value{
a list of three elements: $GEBV is an array of genetic breeding values,
       $phenotypes is the array of phenotypes, and $markerEffects is an array of
       marker effects
}
\description{
Each marker is assigned an effect based on the passed distribution.
Heritability (h2) controls how much GEBV are close to phenotypes (h2=1 means
no distinction, h2=0 means no correlation)
}
\keyword{internal}
