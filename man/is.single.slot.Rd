% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common.R
\name{is.single.slot}
\alias{is.single.slot}
\title{Is a single slot thing?}
\usage{
is.single.slot(x, NULL.is.single = FALSE)
}
\arguments{
\item{x}{the thing to be tested}

\item{NULL.is.single}{should NULL be considered a single slot or not (default: not)}
}
\value{
TRUE if is single slot, FALSE otherwise
}
\description{
Checks if the passed variable is a single slot thing, meaning it contains only a single
value (numeric, character, whatever) and nothing more.
Works with array, vectors, matrix, data.frame...
}
\examples{
is.single.slot(5)   #TRUE
is.single.slot('foobar')   #TRUE
is.single.slot(NULL)       #depends on NULL.is.single
is.single.slot(NA)         #TRUE
is.single.slot(c(1,2,5))   #FALSE
is.single.slot(matrix(0, 10, 5))   #FALSE
is.single.slot(matrix(0, 1, 1))   #TRUE
}
\keyword{internal}
