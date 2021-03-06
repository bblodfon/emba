% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_vector_diff}
\alias{get_vector_diff}
\title{Calculate difference vector with penalty term}
\usage{
get_vector_diff(vec1, vec2, m1 = 1, m2 = 1, penalty = 0)
}
\arguments{
\item{vec1}{numeric vector}

\item{vec2}{numeric vector}

\item{m1}{integer > 0}

\item{m2}{integer > 0}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty (\code{m1,m2} don't matter) and a value of 1 is the strickest possible
penalty. Default value is 0.}
}
\value{
the vector of differences between the two given vectors based on the
formula: \deqn{(vec1 - vec2) * w}, where \eqn{w = (min(m1,m2)/max(m1,m2))^p}
and \eqn{p = penalty}.

See also related \href{https://math.stackexchange.com/questions/3547139/formula-for-weighted-average-difference}{StackOverflow question}.
If \code{vec1} has \code{names}, the returned vector will have the same names
attribute as \code{vec1}.
}
\description{
This function calculates the difference between two given numeric vectors while
adding a penalty term (weight) to account for the number of models/instances that
each vector's values were calculated from. Thus, if the models/instances are
disproportionate and a penalty is included, the difference vector's values will
be changed accordingly to reflect that.
}
