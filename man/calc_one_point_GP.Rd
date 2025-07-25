% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussiantwodRDD.R
\name{calc_one_point_GP}
\alias{calc_one_point_GP}
\title{Gaussian process predictions}
\usage{
calc_one_point_GP(
  sentinels,
  data,
  method = c("new", "aGP", "laGP", "iso"),
  startnum = NULL,
  endnum = NULL,
  prefix = "var."
)
}
\arguments{
\item{sentinels}{Data frame holding the sentinels.}

\item{data}{The data to be analyzed, passed.}

\item{method}{What type of GP to call, passed.}

\item{startnum}{Parameter to Gaussian process estimator.}

\item{endnum}{Parameter to Gaussian process estimator.}

\item{prefix}{Used in Sigma column names.}
}
\value{
Data frame with the sentinels and their GP predictions.  If
variance-covariance available, it will add that in as extra
columns with prefix \code{prefix}.
}
\description{
Innermost function - runs the laGP model to fit gaussian process
and uses it to predict for each (x,y) pair of sentinels.
}
