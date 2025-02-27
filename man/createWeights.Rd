% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createWeights.R
\name{createWeights}
\alias{createWeights}
\title{Create weights for COBRA Predictions}
\usage{
createWeights(train, test, n_train, n_test, nlearners, eps, alpha)
}
\arguments{
\item{train}{A \code{matrix} representing the training data. Rows represent observations,
and columns contain predictions from individual learners for these observations.
In cases where a prediction is unavailable for a specific observation, \code{NA} is used.}

\item{test}{A \code{matrix} representing the test data. Rows represent observations,
and columns contain predictions from individual learners for these observations.
In cases where a prediction is unavailable for a specific observation, \code{NA} is used.}

\item{n_train}{An \code{integer} specifying the number of training observations.}

\item{n_test}{An \code{integer} specifying the number of test observations.}

\item{nlearners}{An \code{integer} representing the number of learners.}

\item{eps}{A \code{numeric} value representing the threshold for proximity between two predictions.}

\item{alpha}{A value that determines the optimal number of learners in the neighborhood (only for alpha optimization).}
}
\description{
The \code{createWeights} function is used to calculate weights for predictions.
}
