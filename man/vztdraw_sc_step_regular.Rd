% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vztdraw_sc_step_regular.R
\name{vztdraw_sc_step_regular}
\alias{vztdraw_sc_step_regular}
\title{Vectorized sampling from zero-truncated NHPPPs with piecewise constant intensities
with same interval lengths}
\usage{
vztdraw_sc_step_regular(
  Lambda_matrix = NULL,
  lambda_matrix = NULL,
  range_t = c(0, 10),
  subinterval = NULL,
  atmost1 = FALSE,
  ...
)
}
\arguments{
\item{Lambda_matrix}{(matrix) integrated intensity rates at the end of each interval}

\item{lambda_matrix}{(matrix) intensity rates, one per interval}

\item{range_t}{(vector, double) \code{t_min} and \code{t_max}}

\item{subinterval}{(vector, double) optional -- the subinterval of \code{range_t} to sample. If \code{NULL}, the whole range_t is used.}

\item{atmost1}{boolean, draw at most 1 event time}

\item{...}{(any) other arguments (ignored  -- used for flexibility in calling from other functions)}
}
\value{
a vector of event times t
if no events realize, it will have 0 length
}
\description{
Simulate a piecewise constant-rate Poisson Point Process over \verb{(t_min, t_max]} (inversion method)
where the intervals have the same length (are "regular").
}
\examples{
x <- vztdraw_sc_step_regular(Lambda_matrix = matrix(1:5, nrow = 1))
}
