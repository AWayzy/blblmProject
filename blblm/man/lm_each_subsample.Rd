% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{lm_each_subsample}
\alias{lm_each_subsample}
\title{compute the estimates}
\usage{
lm_each_subsample(formula, data, n, B, parallel = FALSE, seed = NaN)
}
\arguments{
\item{formula}{The regression formula to be applied}

\item{data}{Dataframe}

\item{n}{Number of rows in the data}

\item{B}{Number of replications to be used in the bootstrap}

\item{parallel}{Determines if the function will utilize multiple cores}

\item{seed}{optional seed for testing}
}
\description{
compute the estimates
}
