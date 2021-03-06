% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/transLocal.R
\name{transLocal}
\alias{transLocal}
\title{Transform the marginals of a multivariate data set to standard normality based
on the logspline density estimator.}
\usage{
transLocal(data, grid = NULL, return.normalizing.constants = FALSE)
}
\arguments{
\item{data}{The data matrix, one row per observation.}

\item{grid}{The grid if it should be transformed as well.}

\item{return.normalizing.constants}{TRUE if the normalizing constants to be
used in the LGDE should be returned as well.}
}
\value{
A list containing the transformed data ($transformed.data), the transformed
grid if provided ($transformed.grid), and the normalizing constants if chosen
($normalizing.constants).
}
\description{
Transform the marginals of a multivariate data set to standard normality based
on the logspline density estimator.
}
\examples{
data <- cbind(rt(100, df = 10), rt(100, df = 10))
grid <- matrix(c(-1, -1, 0, 0, 1, 1), ncol = 2, byrow = TRUE)

transLocal(data, grid, return.normalizing.constants = TRUE)
}

